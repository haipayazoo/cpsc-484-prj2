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
		 double t)
      : _point(point),
	_normal(normal),
	_t(t) {
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
      : SceneObject(diffuse_color, specular_color),
	_center(center),
	_radius(radius) {
      assert(center->is_homogeneous_point());
      assert(radius > 0.0);
    }

    virtual std::shared_ptr<Intersection> intersect(const Vector4& ray_origin,
						    const Vector4& ray_direction) const {
      // See section 4.4.1 of Marschner et al.

      // TODO, this will say objects never intersect.
      return std::shared_ptr<Intersection>(nullptr);
    }
  };

  // Class for a light source.
  class Light {
  private:
    std::shared_ptr<Color> _color;
    double _intensity;

  public:
    Light(std::shared_ptr<Color> color,
	  double intensity)
      : _color(color),
	_intensity(intensity) {
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
    Camera(std::shared_ptr<Vector4> location,
	   std::shared_ptr<Vector4> gaze,
	   std::shared_ptr<Vector4> up,
	   double l, double t, double r, double b,
	   double d)
      : _location(location),
	_gaze(gaze),
	_up(up),
	_l(l),
	_t(t),
	_r(r),
	_b(b),
	_d(d) {
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
      : _ambient_light(ambient_light),
	_background_color(background_color),
	_camera(camera),
	_perspective(perspective) {
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

      // TODO: implement the raytracing algorithm described in section
      // 4.6 of Marschner et al.
      
      return image;
    }

  private:
    // TODO: You will probably want to write some private helper
    // functions to break up the render() function into digestible
    // pieces.
  };
}

// vim: et ts=2 sw=2 :
