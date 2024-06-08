#include <mitsuba/core/bbox.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/render/sensor.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _sensor-perspective:

Perspective pinhole camera with more flexible intrinsic parameter control
(:monosp:`perspective`)
--------------------------------------------------

.. pluginparameters::
 :extra-rows: 7

 * - to_world
   - |transform|
   - Specifies an optional camera-to-world transformation.
     (Default: none (i.e. camera space = world space))
   - |exposed|, |differentiable|, |discontinuous|

 * - fx
   - |Float|
   - Denotes the camera's focal length along the :math:`x`-axis in pixels, e.g.
    K(0, 0) of the camera matrix.
    (Default: :monosp:`width / (2 * dr::tan(0.5 * 45.0))`)
   - |exposed|

 * - fy
   - |Float|
   - Denotes the camera's focal length along the :math:`y`-axis in pixels, e.g.
    K(1, 1) of the camera matrix.
    (Default: :monosp:`fx`)
   - |exposed|

 * - cx
   - |Float|
   - Denotes the camera's principal point along the :math:`x`-axis in pixels,
     to be consistent with graphics community, 0 is defined as the
     corner of the top-left pixel, e.g. K(0, 2) + 0.5 should be used.
    (Default: :monosp:`(width - 1) / 2`)
   - |exposed|

* - cy
   - |Float|
   - Denotes the camera's principal point along the :math:`y`-axis in pixels,
     to be consistent with graphics community, 0 is defined as the
     corner of the top-left pixel, e.g. K(1, 2) + 0.5 of the camera matrix.
    (Default: :monosp:`(height - 1) / 2`)
   - |exposed|

 * - near_clip, far_clip
   - |float|
   - Distance to the near/far clip planes. (Default: :monosp:`near_clip=1e-2`
(i.e. :monosp:`0.01`) and :monosp:`far_clip=1e4` (i.e. :monosp:`10000`))
   - |exposed|

  * - srf
   - |spectrum|
   - Sensor Response Function that defines the :ref:`spectral sensitivity
<explanation_srf_sensor>` of the sensor (Default: :monosp:`none`)
 */

template <typename Float, typename Spectrum>
class IntrinsicCamera final : public ProjectiveCamera<Float, Spectrum> {
public:
    MI_IMPORT_BASE(ProjectiveCamera, m_to_world, m_needs_sample_3, m_film,
                   m_sampler, m_resolution, m_shutter_open, m_shutter_open_time,
                   m_near_clip, m_far_clip, sample_wavelengths)
    MI_IMPORT_TYPES()

    IntrinsicCamera(const Properties &props) : Base(props) {
        ScalarVector2i size = m_film->size();

        m_fx = (ScalarFloat) (size.x() / (2.f * tan(0.3925f)));
        m_fy = m_fx;

        // We record the opengl convention for the principal point, e.g. corner
        // of the corner pixel is 0
        m_cx = (ScalarFloat) ((size.x()) / 2.f);
        m_cy = (ScalarFloat) ((size.y()) / 2.f);

        if (props.has_property("fx")) {
            m_fx = (ScalarFloat) props.get<float>("fx");
        }
        if (props.has_property("fy")) {
            m_fy = (ScalarFloat) props.get<float>("fy");
        }
        if (props.has_property("cx")) {
            m_cx = (ScalarFloat) (props.get<float>("cx"));
        }
        if (props.has_property("cy")) {
            m_cy = (ScalarFloat) (props.get<float>("cy"));
        }

        if (m_to_world.scalar().has_scale())
            Throw("Scale factors in the camera-to-world transformation are not "
                  "allowed!");

        update_camera_transforms();
    }

    void traverse(TraversalCallback *callback) override {
        Base::traverse(callback);
        callback->put_parameter(
            "fx", m_fx, ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_parameter(
            "fy", m_fy, ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_parameter(
            "cx", m_cx, ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_parameter(
            "cy", m_cy, ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_parameter("to_world", *m_to_world.ptr(),
                                ParamFlags::Differentiable |
                                    ParamFlags::Discontinuous);
    }

    void parameters_changed(const std::vector<std::string> &keys) override {
        Base::parameters_changed(keys);
        if (keys.empty() || string::contains(keys, "to_world")) {
            if (m_to_world.scalar().has_scale())
                Throw("Scale factors in the camera-to-world transformation are "
                      "not allowed!");
        }

        update_camera_transforms();
    }

    void update_camera_transforms() {

        Vector2f size = Vector2f(m_film->size());
        Float aspect  = size.x() / size.y();
        Float n       = m_near_clip;
        Float f       = m_far_clip;
        Float recip   = 1.f / (f - n);

        // Basically do the following transform
        // x'=x/z+ox y'=y/z+oy z'=f(z-n)/z(f-n)
        Matrix4f trafo = dr::diag(Vector4f(
            2.f * m_fx / size.x(), 2.f * m_fy / size.y(), f * recip, 0.f));

        trafo(0, 2) = 1.f - 2.f * m_cx / size.x();
        trafo(1, 2) = 1.f - 2.f * m_cy / size.y();
        trafo(2, 3) = -n * f * recip;
        trafo(3, 2) = 1.f;

        /**
         * These do the following (in reverse order):
         *
         * 1. Create transform from camera space to [-1,1]x[-1,1]x[0,1] clip
         *    coordinates (not taking account of the aspect ratio yet)
         *
         * 2+3. Translate and scale to shift the clip coordinates into the
         *    range from zero to one, and take the aspect ratio into account.
         *
         */
        m_camera_to_sample =
            Transform4f::scale(Vector3f(-0.5f, -0.5f * aspect, 1.f)) *
            Transform4f::translate(Vector3f(-1.f, -1.f / aspect, 0.f)) *
            Transform4f(trafo);

        m_sample_to_camera = m_camera_to_sample.inverse();

        // Position differentials on the near plane
        m_dx = m_sample_to_camera * Point3f(1.f / m_resolution.x(), 0.f, 0.f) -
               m_sample_to_camera * Point3f(0.f);
        m_dy = m_sample_to_camera * Point3f(0.f, 1.f / m_resolution.y(), 0.f) -
               m_sample_to_camera * Point3f(0.f);

        /* Precompute some data for importance(). Please
           look at that function for further details. */
        Point3f pmin(m_sample_to_camera * Point3f(0.f, 0.f, 0.f)),
            pmax(m_sample_to_camera * Point3f(1.f, 1.f, 0.f));

        m_image_rect.reset();
        m_image_rect.expand(Point2f(pmin.x(), pmin.y()) / pmin.z());
        m_image_rect.expand(Point2f(pmax.x(), pmax.y()) / pmax.z());
        m_normalization  = 1.f / m_image_rect.volume();
        m_needs_sample_3 = false;

        dr::make_opaque(m_camera_to_sample, m_sample_to_camera, m_fx, m_fy,
                        m_cx, m_cy, m_dx, m_dy, m_image_rect, m_normalization);
    }

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &position_sample,
                                          const Point2f & /*aperture_sample*/,
                                          Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        auto [wavelengths, wav_weight] = sample_wavelengths(
            dr::zeros<SurfaceInteraction3f>(), wavelength_sample, active);
        Ray3f ray;
        ray.time        = time;
        ray.wavelengths = wavelengths;

        // Compute the sample position on the near plane (local camera space).
        Point3f near_p = m_sample_to_camera *
                         Point3f(position_sample.x(), position_sample.y(), 0.f);

        // Convert into a normalized ray direction; adjust the ray interval
        // accordingly.
        Vector3f d = dr::normalize(Vector3f(near_p));

        ray.o = m_to_world.value().translation();
        ray.d = m_to_world.value() * d;

        Float inv_z  = dr::rcp(d.z());
        Float near_t = m_near_clip * inv_z, far_t = m_far_clip * inv_z;
        ray.o += ray.d * near_t;
        ray.maxt = far_t - near_t;

        return { ray, wav_weight };
    }

    std::pair<RayDifferential3f, Spectrum> sample_ray_differential(
        Float time, Float wavelength_sample, const Point2f &position_sample,
        const Point2f & /*aperture_sample*/, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        auto [wavelengths, wav_weight] = sample_wavelengths(
            dr::zeros<SurfaceInteraction3f>(), wavelength_sample, active);
        RayDifferential3f ray;
        ray.time        = time;
        ray.wavelengths = wavelengths;

        // Compute the sample position on the near plane (local camera space).
        Point3f near_p = m_sample_to_camera *
                         Point3f(position_sample.x(), position_sample.y(), 0.f);

        // Convert into a normalized ray direction; adjust the ray interval
        // accordingly.
        Vector3f d = dr::normalize(Vector3f(near_p));

        ray.o = m_to_world.value().translation();
        ray.d = m_to_world.value() * d;

        Float inv_z  = dr::rcp(d.z());
        Float near_t = m_near_clip * inv_z, far_t = m_far_clip * inv_z;
        ray.o += ray.d * near_t;
        ray.maxt = far_t - near_t;

        ray.o_x = ray.o_y = ray.o;

        ray.d_x = m_to_world.value() * dr::normalize(Vector3f(near_p) + m_dx);
        ray.d_y = m_to_world.value() * dr::normalize(Vector3f(near_p) + m_dy);
        ray.has_differentials = true;

        return { ray, wav_weight };
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f & /*sample*/,
                     Mask active) const override {
        // Transform the reference point into the local coordinate system
        Transform4f trafo = m_to_world.value();
        Point3f ref_p     = trafo.inverse().transform_affine(it.p);

        // Check if it is outside of the clip range
        DirectionSample3f ds = dr::zeros<DirectionSample3f>();
        ds.pdf               = 0.f;
        active &= (ref_p.z() >= m_near_clip) && (ref_p.z() <= m_far_clip);
        if (dr::none_or<false>(active))
            return { ds, dr::zeros<Spectrum>() };

        Point3f screen_sample = m_camera_to_sample * ref_p;

        ds.uv = Point2f(screen_sample.x(), screen_sample.y());

        active &= (ds.uv.x() >= 0) && (ds.uv.x() <= 1) && (ds.uv.y() >= 0) &&
                  (ds.uv.y() <= 1);
        if (dr::none_or<false>(active))
            return { ds, dr::zeros<Spectrum>() };

        ds.uv *= m_resolution;

        Vector3f local_d(ref_p);
        Float dist     = dr::norm(local_d);
        Float inv_dist = dr::rcp(dist);
        local_d *= inv_dist;

        ds.p    = trafo.transform_affine(Point3f(0.0f));
        ds.d    = (ds.p - it.p) * inv_dist;
        ds.dist = dist;
        ds.n    = trafo * Vector3f(0.0f, 0.0f, 1.0f);
        ds.pdf  = dr::select(active, Float(1.f), Float(0.f));

        return { ds, Spectrum(importance(local_d) * inv_dist * inv_dist) };
    }

    ScalarBoundingBox3f bbox() const override {
        ScalarPoint3f p = m_to_world.scalar() * ScalarPoint3f(0.f);
        return ScalarBoundingBox3f(p, p);
    }

    /**
     * \brief Compute the directional sensor response function of the camera
     * multiplied with the cosine foreshortening factor associated with the
     * image plane
     *
     * \param d
     *     A normalized direction vector from the aperture position to the
     *     reference point in question (all in local camera space)
     */
    Float importance(const Vector3f &d) const {
        /* How is this derived? Imagine a hypothetical image plane at a
           distance of d=1 away from the pinhole in camera space.

           Then the visible rectangular portion of the plane has the area

              A = (2 * dr::tan(0.5 * xfov in radians))^2 / aspect

           Since we allow crop regions, the actual visible area is
           potentially reduced:

              A' = A * (cropX / filmX) * (cropY / filmY)

           Perspective transformations of such aligned rectangles produce
           an equivalent scaled (but otherwise undistorted) rectangle
           in screen space. This means that a strategy, which uniformly
           generates samples in screen space has an associated area
           density of 1/A' on this rectangle.

           To compute the solid angle density of a sampled point P on
           the rectangle, we can apply the usual measure conversion term:

              d_omega = 1/A' * distance(P, origin)^2 / dr::cos(theta)

           where theta is the angle that the unit direction vector from
           the origin to P makes with the rectangle. Since

              distance(P, origin)^2 = Px^2 + Py^2 + 1

           and

              dr::cos(theta) = 1/sqrt(Px^2 + Py^2 + 1),

           we have

              d_omega = 1 / (A' * cos^3(theta))
        */

        Float ct = Frame3f::cos_theta(d), inv_ct = dr::rcp(ct);

        // Compute the position on the plane at distance 1
        Point2f p(d.x() * inv_ct, d.y() * inv_ct);

        /* Check if the point lies to the front and inside the
           chosen crop rectangle */
        Mask valid = ct > 0 && m_image_rect.contains(p);

        return dr::select(valid, m_normalization * inv_ct * inv_ct * inv_ct,
                          0.f);
    }

    std::string to_string() const override {
        using string::indent;

        std::ostringstream oss;
        oss << "IntrinsicCamera[" << std::endl
            << "  fx = " << m_fx << "," << std::endl
            << "  fy = " << m_fy << "," << std::endl
            << "  cx = " << m_cx << "," << std::endl
            << "  cy = " << m_cy << "," << std::endl
            << "  near_clip = " << m_near_clip << "," << std::endl
            << "  far_clip = " << m_far_clip << "," << std::endl
            << "  film = " << indent(m_film) << "," << std::endl
            << "  sampler = " << indent(m_sampler) << "," << std::endl
            << "  resolution = " << m_resolution << "," << std::endl
            << "  shutter_open = " << m_shutter_open << "," << std::endl
            << "  shutter_open_time = " << m_shutter_open_time << ","
            << std::endl
            << "  to_world = " << indent(m_to_world, 13) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    Transform4f m_camera_to_sample;
    Transform4f m_sample_to_camera;
    BoundingBox2f m_image_rect;
    Float m_normalization;
    Float m_fx, m_fy, m_cx, m_cy;
    Vector3f m_dx, m_dy;
};

MI_IMPLEMENT_CLASS_VARIANT(IntrinsicCamera, ProjectiveCamera)
MI_EXPORT_PLUGIN(IntrinsicCamera, "Intrinsic Camera");
NAMESPACE_END(mitsuba)
