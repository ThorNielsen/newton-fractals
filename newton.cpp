#include <iostream>
#include <complex>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <random>

struct Colour
{
    uint8_t r;
    uint8_t g;
    uint8_t b;
    Colour() : r{0}, g{0}, b{0} {}
    Colour(uint8_t rx, uint8_t gx, uint8_t bx) : r{rx}, g{gx}, b{bx} {}
    Colour(uint32_t c) : r((c&0xFF0000)>>16),
                          g((c&0xFF00)>>8),
                          b(c&0xFF) {}
};

struct Image
{
    std::string name;
    size_t width;
    size_t height;
    std::vector<uint8_t> p;
    Image() = default;
    Image(std::string n, size_t w, size_t h)
    : name{n}, width{w}, height{h}
    {
        p.resize(width*height*3);
    }

    inline Colour pixel(size_t x, size_t y)
    {
        return {p[3*width*y+3*x], p[3*width*y+3*x+1], p[3*width*y+3*x+2]};
    }

    inline void setPixel(size_t x, size_t y, Colour c)
    {
        p[3*width*y+3*x] = c.r;
        p[3*width*y+3*x+1] = c.g;
        p[3*width*y+3*x+2] = c.b;
    }
};

template <typename T>
struct Complex
{
    T real;
    T imag;
    Complex& operator+=(const Complex& o)
    {
        real += o.real;
        imag += o.imag;
        return *this;
    }
    Complex& operator-=(const Complex& o)
    {
        real -= o.real;
        imag -= o.imag;
        return *this;
    }
    Complex& operator*=(const Complex& o)
    {
        T t = real * o.imag + imag * o.real;
        real = real * o.real - imag * o.imag;
        imag += t;
        return *this;
    }
    Complex& operator*=(T f)
    {
        real *= f;
        imag *= f;
        return *this;
    }
    Complex conj() const
    {
        return {real, 0.f - imag};
    }
    T norm() const
    {
        return real * real + imag * imag;
    }
    T arg() const
    {
        return std::atan2(imag, real);
    }
    Complex& operator/=(const Complex& o)
    {
        T denom = T(1) / (o.real * o.real + o.imag * o.imag);
        *this = {real * o.real + imag * o.imag, imag * o.real - real * o.imag};
        real *= denom;
        imag *= denom;
        return *this;
    }
};

template <typename T>
Complex<T> operator+(Complex<T> a, const Complex<T>& b)
{
    return a += b;
}

template <typename T>
Complex<T> operator-(Complex<T> a, const Complex<T>& b)
{
    return a -= b;
}

template <typename T>
Complex<T> operator*(Complex<T> a, const Complex<T>& b)
{
    return {a.real * b.real - a.imag * b.imag,
            a.real * b.imag + a.imag * b.real};
}

template <typename T>
Complex<T> operator*(Complex<T> a, T b)
{
    return a *= b;
}

template <typename T>
Complex<T> operator*(T a, Complex<T> b)
{
    return b *= a;
}

template <typename T>
Complex<T> operator/(Complex<T> a, const Complex<T>& b)
{
    return a /= b;
}

template <typename T>
Complex<T> pow(Complex<T> z, int exp)
{
    Complex<T> ret{T(1), T(0)};
    while (exp > 0)
    {
        if (exp & 1) ret = ret * z;
        z = z * z;
        exp >>= 1;
    }
    return ret;
}

void writePnmHeader(std::ofstream& file, size_t width, size_t height)
{
    file << "P6\n" << width << ' ' << height << "\n255\n";
}

void writeImage(Image img)
{
    std::ofstream file;
    if (img.name[img.name.length() - 4] != '.' ||
        img.name[img.name.length() - 3] != 'p' ||
        img.name[img.name.length() - 2] != 'n' ||
        img.name[img.name.length() - 1] != 'm')
    {
        img.name += ".pnm";
    }
    file.open((img.name).c_str());
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open " +
                                 img.name + ".pnm" +
                                 "for writing.");
    }
    writePnmHeader(file, img.width, img.height);
    if (img.p.size() != 3*img.width*img.height)
    {
        throw std::runtime_error("Image width and/or height is wrong.");
    }
    file.write(reinterpret_cast<char*>(&img.p[0]), img.p.size());
}

// x_{n+1} = x_n - f(x_n)/f'(x_n)

Colour hsv(float h, float s, float v)
{
    float c = s * v;
    h *= 0.954929659f;
    float x = c * (1.f - std::abs(std::fmod(h, 2.f) - 1.f));
    if (h <= 1.f) return Colour(c * 255, x * 255, 0);
    if (h <= 2.f) return Colour(x * 255, c * 255, 0);
    if (h <= 3.f) return Colour(0, c * 255, x * 255);
    if (h <= 4.f) return Colour(0, x * 255, c * 255);
    if (h <= 5.f) return Colour(x * 255, 0, c * 255);
    if (h <= 6.f) return Colour(c * 255, 0, x * 255);
    return Colour(0, 0, 0);
}

template <typename T>
Colour newton(T x, T y, int exp, int maxIt = 256, T maxDist = 1e-12)
{
    Complex<T> z{x, y};
    Complex<T> one{T(1), T(0)};
    Complex<T> p = pow(z, exp);
    Complex<T> pm = pow(z, exp-1);
    int it = 0;
    while ((p-one).norm() > maxDist && ++it < maxIt)
    {
        auto denom = T(exp) * pm;
        z = z - (p-one) / denom;
        p = pow(z, exp);
        pm = pow(z, exp-1);
    }
    return hsv(z.arg() + 3.14159265f, 1.f, 1.f - it / (float)maxIt);// 1. - maxIt / 255.);
}


template <typename T>
Colour newtonAA(T pWidth, T pHeight, int samples,
                const std::vector<T>& positions,
                T x, T y, int exp, int maxIt = 256, T maxDist = 1e-12)
{
    int r = 0, g = 0, b = 0;
    for (int i = 0; i < samples; ++i)
    {
        auto c = newton(x + pWidth*positions[2*i], y+pHeight*positions[2*i+1],
                        exp, maxIt, maxDist);
        r += c.r;
        g += c.g;
        b += c.b;
    }
    return Colour(r / samples, g / samples, b / samples);
}

int main()
{
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist(-.5, .5);
    std::vector<double> positions;
    for (int i = 0; i < 512; ++i)
    {
        positions.push_back(dist(gen));
    }
    float height;
    std::cout << "Complex plane height: ";
    std::cin >> height;
    std::cout << "Image dimensions: ";
    int dimX = 1600, dimY = 900;
    std::cin >> dimX >> dimY;

    int exp;
    std::cout << "Exponent: ";
    std::cin >> exp;
    int samples;
    std::cout << "Anti-aliasing samples: ";
    std::cin >> samples;
    float aspect = dimX * 1.f / dimY;
    float fx = aspect * (height / dimX);
    float px = aspect * height * .5f;
    float fy = height / dimY;
    float py = height * .5f;

    Image img("Newton" + std::to_string(exp)
              + " (" + std::to_string(aspect * height)
              + "x" + std::to_string(height) + ", "
              + std::to_string(samples) + "AA)",
              dimX, dimY);
    for (int x = 0; x < dimX; ++x)
    {
        for (int y = 0; y < dimY; ++y)
        {
            double cx = fx * x - px;
            double cy = fy * y - py;
            img.setPixel(x, y, newtonAA((double)fx, (double)fy, samples,
                                        positions, cx, cy, exp, 256, 1e-16));
        }
    }
    writeImage(img);
}
