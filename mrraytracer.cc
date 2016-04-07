//
// mrraytracer.cc
//
// Demo program that uses the raytracer.hh module.
//
// CPSC 484, CSU Fullerton, Spring 2016, Prof. Kevin Wortman
// Project 2
//
// Name:
//   Kyle Terrien
//   Adam Beck
//   Joe Greene
//
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

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "raytrace.hh"

// Default image dimensions.
const int DEFAULT_WIDTH(640), DEFAULT_HEIGHT(640);

// Hardcoded random number generator seed, for perfect consistency
// between runs.
const int SEED(0xF00DFACE);

// The program can generate one of three hardcoded scenes:
//
// spheres: Trivial scene with only two spheres and one point light
// source.
//
// ballpit: More involved scene with thousands of small balls and
// several point lights.
//
// deluxe: Reserved for expansion; create your own scene!

enum SceneName { SCENE_NAME_SPHERES, SCENE_NAME_BALLPIT, SCENE_NAME_DELUXE };

// Configuration, from command-line arguments.
struct Config {
  SceneName scene_name;
  std::string output_path;
  int width, height;
  bool perspective;
};

// Print command-line usage in the event of user error.
void print_usage() {
  std::cerr << "usage:" << std::endl
	    << "    raytrace_demo [OPTIONS...]" << std::endl
	    << std::endl
	    << "options:" << std::endl
	    << "    --scene SCENE     SCENE must be one of: spheres ballpit deluxe" << std::endl
	    << "    -o OUTPUT_PATH" << std::endl
	    << "    --width W         W must be a positive integer; default is " << DEFAULT_WIDTH << std::endl
	    << "    --height H        H must be a positive integer; default is " << DEFAULT_HEIGHT << std::endl
	    << "    --orthographic    use orthographic projection (default)" << std::endl
	    << "    --perspective     use perspective projection instead of orthographic" << std::endl
	    << std::endl
	    << "Exactly one SCENE and exactly one OUTPUT_PATH must be specified." << std::endl
	    << std::endl;
}

// Convenience function to convert a string to a positive
// integer. Return true on success and false on failure.
bool parse_positive_int(int& result, const std::string& s) {
  try {
    result = std::stoi(s);
    return (result > 0);
  } catch (...) {
    return false;
  }
}

// Parse the command-line arguments to produce an initialized Config
// object.
std::unique_ptr<Config> parse_config(int argc, char** argv) {

  assert(argc > 0);
  assert(argv != nullptr);

  std::vector<std::string> args;
  for (int i = 1; i < argc; ++i) {
    args.push_back(argv[i]);
  }
  
  std::unique_ptr<Config> config(new Config);

  // defaults
  config->width = DEFAULT_WIDTH;
  config->height = DEFAULT_HEIGHT;
  config->perspective = false;

  bool error(false),
    got_scene(false),
    got_output_path(false);

  for (int i = 0; (i < args.size()) && !error; ++i) {
    
    bool last(i == (args.size() - 1));
    
    if (args[i] == "--scene") {
      if (last || got_scene) {
	error = true;
      } else if (args[i+1] == "spheres") {
	config->scene_name = SCENE_NAME_SPHERES;
	i++;
	got_scene = true;
      } else if (args[i+1] == "ballpit") {
	config->scene_name = SCENE_NAME_BALLPIT;
	i++;
	got_scene = true;
      } else if (args[i+1] == "deluxe") {
	config->scene_name = SCENE_NAME_DELUXE;
	i++;
	got_scene = true;
      } else {
	error = true;
      }
    } else if (args[i] == "-o") {
      if (last || got_output_path || args[i+1].empty()) {
	error = true;
      } else {
	config->output_path = args[i+1];
	i++;
	got_output_path = true;
      }
    } else if (args[i] == "--width") {
      if (last || !parse_positive_int(config->width, args[i+1])) {
	error = true;
      } else {
	i++;
      }
    } else if (args[i] == "--height") {
      if (last || !parse_positive_int(config->height, args[i+1])) {
	error = true;
      } else {
	i++;
      }
    } else if (args[i] == "--orthographic") {
      config->perspective = false;
    } else if (args[i] == "--perspective") {
      config->perspective = true;
    } else {
      error = true;
    }
  }

  if (error || !got_scene || !got_output_path) {
    return nullptr;
  } else {
    return config;
  }    
}

int main(int argc, char** argv) {

  auto config(parse_config(argc, argv));
  if (!config) {
    print_usage();
    return 1;
  }

  // Colors to choose from, see
  // http://www.w3schools.com/colors/colors_names.asp
  // for more.
  auto white(raytrace::web_color(0xFFFFFF)),
    near_black(raytrace::web_color(0x202020)),
    pure_red(raytrace::web_color(0xFF0000)),
    pure_green(raytrace::web_color(0x00FF00)),
    pure_blue(raytrace::web_color(0x0000FF)),
    purple(raytrace::web_color(0x800080)),
    orange(raytrace::web_color(0xFFA500)),
    light_yellow(raytrace::web_color(0xFFFFE0)),
    light_blue(raytrace::web_color(0xADD8E6)),
    sky_blue(raytrace::web_color(0x87CEEB)),
    plum(raytrace::web_color(0xDDA0DD)),
    papaya_whip(raytrace::web_color(0xFFEFD5));

  // Reasonable ambient light.
  std::shared_ptr<raytrace::Light> default_ambient_light(new raytrace::Light(light_yellow, 0.25));

  // Origin location.
  auto origin(raytrace::vector4_point(0, 0, 0));

  // Declare a pointer to a scene object.
  std::shared_ptr<raytrace::Scene> scene;

  // Now initialize the scene pointer, based upon the scene specified
  // by command-line arguments.
  switch (config->scene_name) {

  case SCENE_NAME_SPHERES:
    {
      // Two spheres of similar sizes right next to each other.
      
      std::shared_ptr<raytrace::Camera> camera(new raytrace::Camera(origin,
								    raytrace::vector4_translation(0, 0, 1),
								    raytrace::vector4_translation(0, 1, 0),
								    -1, 1,
								    1, -1,
								    2));
  
      scene.reset(new raytrace::Scene(default_ambient_light,
				      near_black,
				      camera,
				      config->perspective));
      
      std::shared_ptr<raytrace::SceneObject> sphere;

      sphere.reset(new raytrace::SceneSphere(plum,
					     white,
					     raytrace::vector4_point(0, 0, 4),
					     0.5));
      scene->add_object(sphere);

      sphere.reset(new raytrace::SceneSphere(papaya_whip,
					     white,
					     raytrace::vector4_point(0.5, 0, 4.5),
					     0.4));
      scene->add_object(sphere);

      std::shared_ptr<raytrace::PointLight> light;
      light.reset(new raytrace::PointLight(white,
					   1.0,
					   raytrace::vector4_point(-2, 1, 0)));
      scene->add_point_light(light);
    }
    break;

  case SCENE_NAME_BALLPIT:
    {
      // We look down on a square playpen of small, brighly colored. spheres.
      
      std::shared_ptr<raytrace::Camera> camera(new raytrace::Camera(raytrace::vector4_point(5, 10, -10),
								    raytrace::vector4_translation(0, -1, 1)->normalized(),
								    raytrace::vector4_translation(0, 1, 0),
								    -1, 1,
								    1, -1,
								    2));
      scene.reset(new raytrace::Scene(default_ambient_light,
				      sky_blue,
				      camera,
				      config->perspective));
      srand(SEED);

      // balls, all coordinates between (0, 0, 0) and (10, 10, 1).
      std::vector<std::shared_ptr<raytrace::Color> > ball_colors;
      ball_colors.push_back(pure_red);
      ball_colors.push_back(pure_green);
      ball_colors.push_back(pure_blue);
      ball_colors.push_back(purple);
      ball_colors.push_back(orange);
      for (int i = 0; i < 8000; ++i) {
	auto color(ball_colors[rand() % ball_colors.size()]);
	double x( (rand() % 1000) / 100.0),
	  y( (rand() % 1000) / 1000.0),
	  z( (rand() % 1000) / 100.0);
	std::shared_ptr<raytrace::SceneObject> object(new raytrace::SceneSphere(color,
										white,
										raytrace::vector4_point(x, y, z),
										0.10));
	scene->add_object(object);
      }

      // lights
      std::shared_ptr<raytrace::PointLight> light;
      light.reset(new raytrace::PointLight(white,
					   1.0,
					   raytrace::vector4_point(-3, 3, 0)));
      scene->add_point_light(light);
      light.reset(new raytrace::PointLight(white,
					   0.8,
					   raytrace::vector4_point(1, 3, 0)));
      scene->add_point_light(light);
      light.reset(new raytrace::PointLight(white,
					   0.4,
					   origin));
      scene->add_point_light(light);
    }
    break;
    
  default:
    std::cerr << "ERROR: sorry, that scene is not supported" << std::endl;
    return 1;
  }

  // Check that the scene pointer really did get initialized.
  assert(scene != nullptr);

  // Raytrace!
  auto image(scene->render(config->width, config->height));
  if (!image) {
    std::cerr << "ERROR: rendering error" << std::endl;
    return 1;
  }

  // Write the image to disk.
  if (!image->write_ppm(config->output_path)) {
    std::cerr << "ERROR: could not write " << config->output_path << std::endl;
    return 1;
  }

  // Success.
  return 0;
}

// vim: et ts=2 sw=2 :
