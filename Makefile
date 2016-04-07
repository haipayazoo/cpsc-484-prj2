CFLAGS := -Wall -std=c++11 -Wextra -Wpedantic
CC := clang++

mrraytracer: gmath.hh raytrace.hh mrraytracer.cc
	$(CC) $(CFLAGS) mrraytracer.cc -o mrraytracer

clean:
	-rm -f mrraytracer spheres_o.ppm spheres_p.ppm ballpit_o.ppm ballpit_p.ppm

all: mrraytracer

test: spheres_o.ppm

images: spheres_o.ppm spheres_p.ppm ballpit_o.ppm ballpit_p.ppm

spheres_o.ppm: mrraytracer
	./mrraytracer --scene spheres -o spheres_o.ppm

spheres_p.ppm: mrraytracer
	./mrraytracer --scene spheres -o spheres_p.ppm --perspective

ballpit_o.ppm: mrraytracer
	./mrraytracer --scene ballpit -o ballpit_o.ppm

ballpit_p.ppm: mrraytracer
	./mrraytracer --scene ballpit -o ballpit_p.ppm --perspective

.PHONY: all clean
