#ifndef __CUSTOMCONTAINERS_HPP__
#define __CUSTOMCONTAINERS_HPP__

#include <cstdlib>

namespace CustomContainer {
    class Direction3D {
        private:
            double dx, dy, dz;

        public:    
            Direction3D(double dx0 = 0., double dy0 = 0., double dz0 = 0.): dx(dx0), dy(dy0), dz(dz0) {}

            void setValue(double dx0, double dy0, double dz0){
                this->dx = dx0;
                this->dy = dy0;
                this->dz = dz0;
            }

            double GetDX() { return this->dx; }

            double GetDY() { return this->dy; }

            double GetDZ() { return this->dz; }

            ClassDef(Direction3D, 1);
        };

    class Position4D {
        private:
            double x, y, z, t;

        public:
            Position4D(double x0 = 0., double y0 = 0., double z0 = 0., double t0 = 0.): x(x0), y(y0), z(z0), t(t0) {}

            void setValue(double x0, double y0, double z0, double t0) {
                this->x = x0;
                this->y = y0;
                this->z = z0;
                this->t = t0;
            }

            double GetX() { return this->x; }

            double GetY() { return this->y; }
            
            double GetZ() { return this->z; }

            double GetT() { return this->t; }

            ClassDef(Position4D, 1);
    };

}

#endif