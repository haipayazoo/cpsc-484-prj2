//
// raytrace.hh
//
// Basic raytracer module.
//
// CPSC 484, CSU Fullerton, Spring 2016, Prof. Kevin Wortman
// Project 2
//
// Name:
//   Kyle Terrien
//   Adam Beck
//   Joe Greene
///
// In case it ever matters, this file is hereby placed under the MIT
// License:
//
// Copyright (c) 2016, Kevin Wortman
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#pragma once

#include <cassert>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

#include "gmath.hh"

namespace raytrace {

  // Type aliases for color, vector, and matrix. A color is always a
  // 3-dimensional R, G, B vector; we always use 4-dimensional vectors
  // and matrices for 3D coordinates with homogeneous coordinates.
  
  typedef gmath::Vector<double, 3> Color;
  
  typedef gmath::Vector<double, 4> Vector4;

  typedef gmath::Matrix<double, 4, 4> Matrix4x4;

  // Convenience functions to create vectors.
  
  std::shared_ptr<Vector4> vector4(double x, double y, double z, double w) {
    std::shared_ptr<Vector4> v(new Vector4);
    (*v)[0] = x;
    (*v)[1] = y;
    (*v)[2] = z;
    (*v)[3] = w;
    return v;
  }

  std::shared_ptr<Vector4> vector4_point(double x, double y, double z) {
    return vector4(x, y, z, 1.0);
  }

  std::shared_ptr<Vector4> vector4_translation(double x, double y, double z) {
    return vector4(x, y, z, 0.0);
  }

  // Test whether a scalar represents an R, G, or B intensity in the
  // range [0, 1].
  
  bool is_color_intensity(double x) {
    return ((x >= 0.0) && (x <= 1.0));
  }

  // Test whether a 3-vector represents a valid R, G, B color.
  
  bool is_color(const Color& c) {
    return (is_color_intensity(c[0]) &&
      is_color_intensity(c[1]) &&
      is_color_intensity(c[2]));
  }

  // Convenience function to convert a 24-bit hexadecimal web color,
  // as used in HTML, to one of our Color objects.
  std::shared_ptr<Color> web_color(uint_fast32_t hex) {
    assert(hex <= 0xFFFFFF);
    std::shared_ptr<Color> color(new Color);
    (*color)[0] = (hex >> 16) / 255.0;
    (*color)[1] = ((hex >> 8) & 0xFF) / 255.0;
    (*color)[2] = (hex & 0xFF) / 255.0;
    assert(is_color(*color));
    return color;
  }

  // Class that represents an intersection between a viewing ray and a
  // scene object, complete with a point of intersection, surface
  // normal vector, and time parameter t.
  class Intersection {
  private:
    std::shared_ptr<Vector4> _point, _normal;
    double _t;

  public:
    Intersection(std::shared_ptr<Vector4> point,
        std::shared_ptr<Vector4> normal,
        double t
        ) : _point(point), _normal(normal), _t(t) {
      assert(point->is_homogeneous_point());
      assert(normal->is_homogeneous_translation());
      assert(t >= 0.0);
    }

    const Vector4& point() const { return *_point; }
    const Vector4& normal() const { return *_normal; }
    double t() const { return _t; }
  };

  // Abstract class for a scene object. In a production raytracer we'd
  // have many subclasses for spheres, planes, triangles, meshes,
  // etc. For now we will only have one subclass representing a
  // sphere.
  class SceneObject {
  private:
    std::shared_ptr<Color> _diffuse_color, _specular_color;

  public:
    SceneObject(std::shared_ptr<Color> diffuse_color,
    std::shared_ptr<Color> specular_color)
      : _diffuse_color(diffuse_color),
  _specular_color(specular_color) {
      assert(is_color(*diffuse_color));
      assert(is_color(*specular_color));
    }

    const Color& diffuse_color() const { return *_diffuse_color; }
    const Color& specular_color() const { return *_specular_color; }

    // Abstract virtual function for intersection testing. Given a
    // viewing ray defined by an origin and direction, return a
    // pointer to an Intersection object representing where the ray
    // intersects with the object. If they never intersect, return
    // nullptr.
    virtual std::shared_ptr<Intersection> intersect(const Vector4& ray_origin,
                const Vector4& ray_direction) const = 0;
  };

  // Concrete subclass for a sphere.
  class SceneSphere : public SceneObject {
  private:
    std::shared_ptr<Vector4> _center;
    double _radius;

  public:
    SceneSphere(std::shared_ptr<Color> diffuse_color,
    std::shared_ptr<Color> specular_color,
    std::shared_ptr<Vector4> center,
    double radius)
      : SceneObject(diffuse_color, specular_color), _center(center), _radius(radius) {
      assert(center->is_homogeneous_point());
      assert(radius > 0.0);
    }

    virtual std::shared_ptr<Intersection> intersect(const Vector4& ray_origin,
                const Vector4& ray_direction) const {
      // See section 4.4.1 of Marschner et al.
      // reference: 
      // http://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

      // TODO: review below code (unsure if it actually is complete/valid)

      // Geometric way of computing intersection between a ray and sphere

      double t_m = ray_origin * ray_direction * -1;
      double l_m2 = ray_origin * ray_direction - (ray_origin * ray_direction) * (ray_origin * ray_direction);
      
      // from book (we think; double-check later); we'd have to make sure that l_m2 <= 1 too
      if(l_m2 < 0)
      {
        return nullptr;
      }
      double delta_t = sqrt(1 - l_m2);
      double t0 = -t_m + delta_t;
      double t1 = -t_m - delta_t;

      // (xor) swap (mimics work done in scratchapixel.com tutorial)
      if(t0 > t1) {
        double temp = t0;
        t0 = t1;
        t1 = temp;
      }

      // Intersection (use positive t, as noted in the scratchapixel.com tutorial)
      // NOTE: Earlier, on l_m2 < 0 check, we already verified that an intersection exists.
      //       Thus, no need to worry about a "bad" t here.

      /* 

          QUESTION: Is ray_origin the point and ray_direction the normal?
            Idea 1:
              // time
              time = (t0 < 0) ? t1 : t0;
              // hit_point
              hit = p(t) = ray_origin + ray_direction * t;
              // hit_normal
              hit_normal = hit / hit->magnitude();

              // Intersection
              return std::shared_ptr<Intersection>(new Intersection(hit, hit.normal(), time));
            Idea 2:
              // time
              time = (t0 < 0) ? t1 : t0;
              // hit point: using p + time * d
              hit_point = ray_origin + (ray_direction * time);
              hit_normal = (hit_point - _center) * 2;

      */

      // POSSIBLE TODO: Create #define make_vector_ptr(i) (v*i)

      return std::shared_ptr<Intersection>(
          new Intersection(ray_origin*1,         // lazy way to return a ptr type of ray_origin
                           ray_direction*1,      // "                              " ray_direction
                           (t0 < 0) ? t1 : t0)); // choose less-negative t-value
    }
  };

  // Class for a light source.
  class Light {
  private:
    std::shared_ptr<Color> _color;
    double _intensity;

  public:
    Light(std::shared_ptr<Color> color, double intensity)
      : _color(color), _intensity(intensity) {
      assert(is_color(*color));
      assert(intensity > 0.0);
    }

    const Color& color() const { return *_color; }
    double intensity() const { return _intensity; }
  };

  // A PointLight is a subclass of Light that adds a location in
  // space.
  class PointLight : public Light {
  private:
    std::shared_ptr<Vector4> _location;

  public:
    PointLight(std::shared_ptr<Color> color,
         double intensity,
         std::shared_ptr<Vector4> location)
      : Light(color, intensity),
  _location(location) {
      assert(location->is_homogeneous_point());
    }
      
    const Vector4& location() const { return *_location; }
  };

  // A camera, defined by a location; gaze vector; up vector; viewing
  // plane bounds l, t, r, b; and viewing plane distance d.
  class Camera {
  private:
    std::shared_ptr<Vector4> _location, _gaze, _up;
    double _l, _t, _r, _b, _d;

  public:
    Camera(std::shared_ptr<Vector4> location, std::shared_ptr<Vector4> gaze,
     std::shared_ptr<Vector4> up, double l, double t, double r, double b, double d)
      : _location(location), _gaze(gaze), _up(up), _l(l), _t(t), _r(r), _b(b), _d(d) {
      assert(location->is_homogeneous_point());
      assert(gaze->is_homogeneous_translation());
      assert(up->is_homogeneous_translation());
      assert((l < 0.0) && (0.0 < r));
      assert((b < 0.0) && (0.0 < t));
      assert(d > 0.0);
    }

    const Vector4& location() const { return *_location; }
    const Vector4& gaze() const { return *_gaze; }
    const Vector4& up() const { return *_up; }
    double l() const { return _l; }
    double t() const { return _t; }
    double r() const { return _r; }
    double b() const { return _b; }
    double d() const { return _d; }
  };

  // A raster image, i.e. a rectangular grid of Color objects.
  class Image {
  private:
    std::vector<std::vector<Color> > _pixels;

  public:
    // Initialize the image with the given width and height, and every
    // pixel initialized to fill.
    Image(int width, int height, const Color& fill) {
      assert(width > 0);
      assert(height > 0);
      _pixels.assign(height, std::vector<Color>(width, fill));
    }

    int width() const { return _pixels[0].size(); }
    int height() const { return _pixels.size(); }

    // Determine whether a given int is a valid x/y coordinate.
    bool is_x_coordinate(int x) const {
      return ((x >= 0) && (x < width()));
    }
    bool is_y_coordinate(int y) const {
      return ((y >= 0) && (y < height()));
    }

    // Determine whether (x, y) ints represent valid coordinates.
    bool is_coordinate(int x, int y) const {
      return (is_x_coordinate(x) && is_y_coordinate(y));
    }

    // Get or set a single pixel.
    const Color& pixel(int x, int y) const {
      assert(is_coordinate(x, y));
      return _pixels[y][x];
    }
    void set_pixel(int x, int y, const Color& color) {
      assert(is_coordinate(x, y));
      assert(is_color(color));
      _pixels[y][x] = color;
    }

    // Write the image to a file in the PPM file format.
    // https://en.wikipedia.org/wiki/Netpbm_format
    //
    // Return true on success or false in the case of an I/O error.
    bool write_ppm(const std::string& path) const {
      std::ofstream f(path);
      if (!f)
        return false;

      f << "P3" << std::endl
        << width() << ' ' << height() << std::endl
        << "255" << std::endl;

      for (int y = height()-1; y >= 0; --y) {
        for (int x = 0; x < width(); ++x) {
          const Color& c(pixel(x, y));
          if (x > 0)
            f << ' ';
            f << discretize(c[0]) << ' '
              << discretize(c[1]) << ' '
              << discretize(c[2]);
        }
        f << std::endl;
      }

      bool success(f);
      f.close();
      
      return success;
    }

  private:
    // Convert a scalar color intensity in the range [0, 1] to a byte
    // value in the range [0, 255].
    int discretize(double intensity) const {
      assert(is_color_intensity(intensity));
      int x = static_cast<int>(round(intensity * 255.0));
      if (x < 0)
        x = 0;
      else if (x > 255)
        x = 255;
      return x;
    }
  };

  // Class for an entire scene, tying together all the other classes
  // in this module.
  class Scene {
  private:
    // Ambient light source, to prevent objects that are blocked from
    // point light sources from being entirely black.
    std::shared_ptr<Light> _ambient_light;

    // Background color for pixels that do not correspond to any scene
    // object.
    std::shared_ptr<Color> _background_color;

    // The camera.
    std::shared_ptr<Camera> _camera;

    // When true, use the perspective projection; when false, use
    // orthographic projection.
    bool _perspective;

    // Vector of all scene objects.
    std::vector<std::shared_ptr<SceneObject>> _objects;

    // Vector of all point lights.
    std::vector<std::shared_ptr<PointLight>> _point_lights;

  public:
    // Initialize a scene, initially with no objects and no point
    // lights.
    Scene(std::shared_ptr<Light> ambient_light,
        std::shared_ptr<Color> background_color,
        std::shared_ptr<Camera> camera,
        bool perspective)
      : _ambient_light(ambient_light), _background_color(background_color),
      _camera(camera), _perspective(perspective) {
      assert(is_color(*background_color));
    }

    // Add an object/light.
    void add_object(std::shared_ptr<SceneObject> object) { _objects.push_back(object); }
    void add_point_light(std::shared_ptr<PointLight> light) { _point_lights.push_back(light); }

    // Render the scene into an image of the given width and height.
    //
    // This is the centerpiece of the module, and is responsible for
    // executing the core raytracing algorithm.
    std::shared_ptr<Image> render(int width, int height) const {
      assert(width > 0);
      assert(height > 0);
      
      std::shared_ptr<Image> image(new Image(width, height, *_background_color));

      // pixel coordinate positions
      int i, j;

      // viewing ray pointers
      std::shared_ptr<Vector4> ray_origin, ray_direction;

      // for each pixel
      for (j = 0; j < height; ++j) {
        for (i = 0; i < width; ++i) {
          for (auto scene_obj : _objects) {
            // compute viewing ray (initializes ray_origin and ray_direction)
            compute_viewing_ray(ray_origin, ray_direction);

            // compute intersection point
            std::shared_ptr<Intersection> hit_point = scene_obj->intersect(ray_origin, ray_direction);

            // if an intersection exists between the viewing ray and scene object
            if(hit_point != nullptr) {
              // TODO: Review below work
              std::shared_ptr<Vector4> surface_normal = hit_point->normal();

              // evaluate shading model and set pixel to that color // page 82
              // we have diffuse_color and light intensity
              //image->set_pixel(i,j, light_intensity());
            }
            else
            {
              // set pixel color to background color (no hit)
              //image->set_pixel(i, j, _background_color);
            }
          }
        }
      }
      // TODO: implement the raytracing algorithm described in section
      // 4.6 of Marschner et al.
      
      return image;
    }

  private:
    // computes viewing ray
    void compute_viewing_ray(std::shared_ptr<Vector4>& ray_origin, 
                             std::shared_ptr<Vector4>& ray_direction){
      // what are these again?
      double u, v;

      // what are these again? unit vectors of the camera or something like that
      std::shared_ptr<Vector4> vec_u, vec_v, vec_w;

      // compute u and v
      u = _camera->l() + (_camera->r() - _camera->l()) * (i + 0.5) / width;
      v = _camera->b() + (_camera->t() - _camera->b()) * (j + 0.5) / height;

      // if we are using the perspective transform, compute the viewing ray based off of 
      //     perspective transform
      if(_perspective) {
        // TODO: somehow clean up below line
        ray_direction = *(*(*vec_w * -_camera->d()) + *(*vec_u * u)) + *(*vec_v * v);

        // saying "vector4 * 1" is a lazy way to make a ptr_type of the vector4
        // POSSIBLE TODO: create macro #define make_vector_ptr(i) (v*i)
        ray_origin = _camera->location() * 1;
      }
      else {
        ray_direction = -(*vec_w);
        // TODO: somehow clean up below line
        ray_origin = *(_camera->location() + *(*vec_u * u)) + *vec_v * v;
      }
    }
    /*
    // unsure about return type
    void evaluateShading(){
      
      // set pixel to correct color
    
    }*/
    // TODO: You will probably want to write some private helper
    // functions to break up the render() function into digestible
    // pieces.
  };
}

// vim: et ts=2 sw=2 :
