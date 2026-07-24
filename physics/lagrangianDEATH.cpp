/******************************************************************************************
*
*  COMPUTATIONAL PHYSICS ENGINE
*
*  Demonstration of a small object-oriented spring-mass simulation.
*
*  Features
*  --------
*  • Vector algebra
*  • Particles
*  • Springs
*  • Gravity
*  • Damping
*  • RK4 Integrator
*  • Energy computation
*  • CSV Export
*
******************************************************************************************/

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

constexpr double PI = 3.14159265358979323846;

class Vector2
{
public:

    double x;
    double y;

    Vector2()
        : x(0.0), y(0.0)
    {
    }

    Vector2(double X, double Y)
        : x(X), y(Y)
    {
    }

    Vector2 operator+(const Vector2& rhs) const
    {
        return Vector2(x + rhs.x, y + rhs.y);
    }

    Vector2 operator-(const Vector2& rhs) const
    {
        return Vector2(x - rhs.x, y - rhs.y);
    }

    Vector2 operator*(double s) const
    {
        return Vector2(x * s, y * s);
    }

    Vector2 operator/(double s) const
    {
        return Vector2(x / s, y / s);
    }

    Vector2& operator+=(const Vector2& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    double magnitude() const
    {
        return std::sqrt(x*x + y*y);
    }

    Vector2 normalized() const
    {
        double m = magnitude();

        if(m < 1e-10)
            return Vector2();

        return *this / m;
    }

    static double dot(
        const Vector2& a,
        const Vector2& b)
    {
        return a.x*b.x + a.y*b.y;
    }
};

class Particle
{
public:

    double mass;

    Vector2 position;

    Vector2 velocity;

    Vector2 accumulatedForce;

    Particle(
        double m,
        Vector2 pos
    )
        :
        mass(m),
        position(pos),
        velocity(),
        accumulatedForce()
    {
    }

    void clearForce()
    {
        accumulatedForce = Vector2();
    }

    void applyForce(const Vector2& F)
    {
        accumulatedForce += F;
    }

    Vector2 acceleration() const
    {
        return accumulatedForce / mass;
    }

    double kineticEnergy() const
    {
        return 0.5 *
               mass *
               Vector2::dot(
                    velocity,
                    velocity
               );
    }
};

class Spring
{
public:

    Particle* a;

    Particle* b;

    double stiffness;

    double restLength;

    Spring(
        Particle* A,
        Particle* B,
        double K,
        double L
    )
        :
        a(A),
        b(B),
        stiffness(K),
        restLength(L)
    {
    }

    void apply()
    {
        Vector2 displacement =
            b->position - a->position;

        double length =
            displacement.magnitude();

        if(length < 1e-8)
            return;

        Vector2 direction =
            displacement.normalized();

        double extension =
            length - restLength;

        Vector2 force =
            direction *
            (stiffness * extension);

        a->applyForce(force);

        b->applyForce(force * (-1));
    }

    double potentialEnergy() const
    {
        Vector2 displacement =
            b->position - a->position;

        double extension =
            displacement.magnitude()
            - restLength;

        return
            0.5 *
            stiffness *
            extension *
            extension;
    }
};

class PhysicsWorld
{
private:

    std::vector<
        std::unique_ptr<Particle>
    > particles;

    std::vector<
        std::unique_ptr<Spring>
    > springs;

    Vector2 gravity;

    double damping;

public:

    PhysicsWorld()
        :
        gravity(0.0,-9.81),
        damping(0.25)
    {
    }

    Particle* createParticle(
        double m,
        Vector2 p
    )
    {
        particles.push_back(
            std::make_unique<
                Particle
            >(m,p)
        );

        return particles.back().get();
    }

    void createSpring(
        Particle* a,
        Particle* b,
        double k,
        double L
    )
    {
        springs.push_back(
            std::make_unique<
                Spring
            >(a,b,k,L)
        );
    }

    void computeForces()
    {
        for(auto& p : particles)
        {
            p->clearForce();

            p->applyForce(
                gravity *
                p->mass
            );

            p->applyForce(
                p->velocity *
                (-damping)
            );
        }

        for(auto& s : springs)
        {
            s->apply();
        }
    }

    void integrateRK4(double dt)
    {
        computeForces();

        for(auto& p : particles)
        {
            Vector2 a =
                p->acceleration();

            Vector2 k1v = a;
            Vector2 k1x = p->velocity;

            Vector2 k2v = a;
            Vector2 k2x =
                p->velocity +
                k1v*(dt*0.5);

            Vector2 k3v = a;
            Vector2 k3x =
                p->velocity +
                k2v*(dt*0.5);

            Vector2 k4v = a;
            Vector2 k4x =
                p->velocity +
                k3v*dt;

            p->position +=
                (k1x +
                 k2x*2.0 +
                 k3x*2.0 +
                 k4x)
                 *
                 (dt/6.0);

            p->velocity +=
                (k1v +
                 k2v*2.0 +
                 k3v*2.0 +
                 k4v)
                 *
                 (dt/6.0);
        }
    }

    double totalEnergy() const
    {
        double E = 0.0;

        for(const auto& p : particles)
        {
            E +=
                p->kineticEnergy();
        }

        for(const auto& s : springs)
        {
            E +=
                s->potentialEnergy();
        }

        return E;
    }

    void exportCSV(
        const std::string& filename,
        int steps,
        double dt)
    {
        std::ofstream file(filename);

        file
        << "time,"
        << "x,"
        << "y,"
        << "energy\n";

        for(int i=0;i<steps;i++)
        {
            integrateRK4(dt);

            file
            << std::fixed
            << std::setprecision(6);

            file
            << i*dt
            << ","
            << particles[0]->position.x
            << ","
            << particles[0]->position.y
            << ","
            << totalEnergy()
            << "\n";
        }

        file.close();
    }
};

int main()
{
    PhysicsWorld world;

    Particle* p1 =
        world.createParticle(
            1.0,
            {0.0,0.0}
        );

    Particle* p2 =
        world.createParticle(
            1.0,
            {2.0,0.0}
        );

    world.createSpring(
        p1,
        p2,
        100.0,
        1.5
    );

    world.exportCSV(
        "simulation.csv",
        5000,
        0.001
    );

    std::cout
        << "Simulation Complete\n";

    return 0;
}
