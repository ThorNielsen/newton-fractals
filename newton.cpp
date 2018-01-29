#include <iostream>
#include <complex>
#include <stdexcept>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <thread>

#include "SFML/Window.hpp"
#include "SFML/Graphics.hpp"

struct Colour
{
    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t a;
    Colour() : r{0}, g{0}, b{0}, a{255} {}
    Colour(uint8_t rx, uint8_t gx, uint8_t bx) : r{rx}, g{gx}, b{bx}, a{255} {}
    Colour(uint32_t c) : r((c&0xFF0000)>>16),
                          g((c&0xFF00)>>8),
                          b(c&0xFF),
                          a{255} {}
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
        p.resize(width*height*4);
    }

    void resize(size_t w, size_t h)
    {
        p.resize(w*h*4);
        width = w;
        height = h;
    }

    inline Colour pixel(size_t x, size_t y)
    {
        return {p[4*width*y+4*x], p[4*width*y+4*x+1], p[4*width*y+4*x+2]};
    }

    inline void setPixel(size_t x, size_t y, Colour c)
    {
        p[4*width*y+4*x] = c.r;
        p[4*width*y+4*x+1] = c.g;
        p[4*width*y+4*x+2] = c.b;
        p[4*width*y+4*x+3] = 255;
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

// x_{n+1} = x_n - f(x_n)/f'(x_n)
template <typename T>
Colour newton(Complex<T> chaos, T x, T y, int exp, int maxIt = 256,
              T maxDist = 1e-12)
{
    Complex<T> z{x, y};
    Complex<T> one{T(1), T(0)};
    Complex<T> p = pow(z, exp);
    Complex<T> pm = pow(z, exp-1);
    int it = 0;
    auto pmn = (p-one).norm();
    // Gives an interesting divergent effect:
    // while (pmn > maxDist && ++it < maxIt && pmn <= std::numeric_limits<T>::max())
    while (pmn > maxDist && ++it < maxIt
           && p.real <= std::numeric_limits<T>::max()
           && p.imag <= std::numeric_limits<T>::max())
    {
        auto denom = T(exp) * pm;
        z = z - chaos * (p-one) / denom;
        p = pow(z, exp);
        pm = pow(z, exp-1);
        pmn = (p-one).norm();
    }
    return hsv(z.arg() + 3.14159265f, 1.f, 1.f - 0.05f * it / (float)maxIt);
}


template <typename T>
Colour newtonAA(T pWidth, T pHeight, int samples,
                const std::vector<T>& positions, Complex<T> chaos,
                T x, T y, int exp, int maxIt = 256, T maxDist = 1e-12)
{
    int r = 0, g = 0, b = 0;
    for (int i = 0; i < samples; ++i)
    {
        auto c = newton(chaos, x + pWidth*positions[2*i],
                        y+pHeight*positions[2*i+1], exp, maxIt, maxDist);
        r += c.r;
        g += c.g;
        b += c.b;
    }
    return Colour(r / samples, g / samples, b / samples);
}

struct FractalInfo
{
    Complex<double> chaos;
    double width;
    double height;
    int exponent;
    int samples;
    int maxIterations;
    int pixelWidth;
    int pixelHeight;
    void recalculateWidth()
    {
        double aspect = static_cast<double>(pixelWidth)
                        / static_cast<double>(pixelHeight);
        width = aspect * height;
    }
    void recalculateHeight()
    {
        double aspect = static_cast<double>(pixelHeight)
                        / static_cast<double>(pixelWidth);
        height = aspect * width;
    }
};

void calculateNewtonChunk(int xBegin, int xEnd, int yBegin, int yEnd,
                          FractalInfo info, Image& img,
                          const std::vector<double>& aaPos,
                          int* done)
{
    double fx = info.width / static_cast<double>(info.pixelWidth);
    double px = info.width * 0.5;
    double fy = info.height / static_cast<double>(info.pixelHeight);
    double py = info.height * 0.5;
    for (int x = xBegin; x < xEnd; ++x)
    {
        for (int y = yBegin; y < yEnd; ++y)
        {
            double cx = fx * x - px;
            double cy = fy * y - py;
            img.setPixel(x, y, newtonAA(fx, fy, info.samples, aaPos, info.chaos,
                                        cx, cy, info.exponent,
                                        info.maxIterations, 1e-16));
        }
    }
    std::cout << "Calculated chunk [" << xBegin << ", " << xEnd << ")×["
              << yBegin << ", " << yEnd << ")" << std::endl;
    *done = 1;
}

FractalInfo getInfoFromStdin()
{
    FractalInfo info;
    std::cout << "Complex plane height: ";
    std::cin >> info.height;
    std::cout << "Image dimensions: ";
    std::cin >> info.pixelWidth >> info.pixelHeight;
    std::cout << "Exponent: ";
    std::cin >> info.exponent;
    std::cout << "Anti-aliasing samples: ";
    std::cin >> info.samples;
    std::cout << "Max iterations: ";
    std::cin >> info.maxIterations;
    std::cout << "Chaos: ";
    std::cin >> info.chaos.real >> info.chaos.imag;
    info.recalculateWidth();
    return info;
}

FractalInfo defaultInfo()
{
    FractalInfo info;
    info.height = 6.;
    info.pixelWidth = 800;
    info.pixelHeight = 600;
    info.exponent = 6;
    info.samples = 51;
    info.maxIterations = 64;
    info.chaos.real = 22.;
    info.chaos.imag = 19.;
    info.recalculateWidth();
    return info;
}

std::string title(FractalInfo info)
{
    return "Newton fractal " + std::to_string(info.exponent)
           + " (" + std::to_string(info.width)
           + "x" + std::to_string(info.height) + ", "
           + std::to_string(info.samples) + "AA, chaos:"
           + std::to_string(info.chaos.real)
           + (info.chaos.imag >= 0. ? "+" : "")
           + std::to_string(info.chaos.imag)
           + "i, "
           + std::to_string(info.maxIterations)
           + " iterations)";
}

std::vector<double> getAAPos(size_t samples)
{
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist(-.5, .5);
    std::vector<double> positions;
    for (int i = 0; i < 2 * samples; ++i)
    {
        positions.push_back(dist(gen));
    }
    return positions;
}

struct ImageRenderer
{
    Image img;
    FractalInfo info;
    ~ImageRenderer()
    {
        for (auto& thread : threads)
        {
            thread.join();
        }
    }

    bool isReady() const
    {
        size_t i = 0;
        for (auto d : done)
        {
            if (!d) return false;
        }
        return true;
    }

    void calculate()
    {
        if (!isReady()) return;
        done.clear();
        for (auto& thread : threads)
        {
            thread.join();
        }
        threads.clear();

        std::cout << "Resizing image to " << info.pixelWidth << "x" << info.pixelHeight << std::endl;
        img.resize(info.pixelWidth, info.pixelHeight);

        int yBegin = 0;
        auto maxThreads = std::thread::hardware_concurrency();
        done.resize(maxThreads);

        int chunkHeight = info.pixelHeight / maxThreads;
        auto positions = getAAPos(info.samples);
        for (unsigned i = 0; i < maxThreads; ++i)
        {
            threads.emplace_back(calculateNewtonChunk, 0, info.pixelWidth,
                                 yBegin, yBegin + chunkHeight, info, std::ref(img),
                                 positions, &done[i]);
            yBegin += chunkHeight;
        }

        int dummy = 0;
        calculateNewtonChunk(0, info.pixelWidth, yBegin, info.pixelHeight,
                             info, img, positions, &dummy);
        std::cout << dummy << std::endl;
    }

private:
    std::vector<int> done;
    std::vector<std::thread> threads;
};

int main(int argc, char* argv[])
{
    ImageRenderer ir;
    ir.info = defaultInfo();
    if (argc > 1)
    {
        ir.info = getInfoFromStdin();
    }

    sf::RenderWindow wnd(sf::VideoMode(ir.info.pixelWidth, ir.info.pixelHeight), "Bob");

    ir.calculate();

    decltype(wnd.getSize()) lastSize = {0, 0};

    while (wnd.isOpen())
    {
        sf::Event evt;
        while (wnd.pollEvent(evt))
        {
            if (evt.type == sf::Event::Closed)
            {
                wnd.close();
            }
        }
        wnd.clear(sf::Color::Blue);
        auto sz = wnd.getSize();
        wnd.setView(sf::View(sf::FloatRect(0, 0, sz.x, sz.y)));
        if (lastSize != sz && ir.isReady())
        {
            lastSize = sz;
            ir.info.pixelWidth = sz.x;
            ir.info.pixelHeight = sz.y;
            ir.info.recalculateWidth();
            ir.calculate();
        }
        sf::Texture texture;
        texture.create(ir.info.pixelWidth, ir.info.pixelHeight);
        texture.update(ir.img.p.data());
        sf::Sprite sprite;
        sprite.setTexture(texture);

        wnd.draw(sprite);
        wnd.display();
    }
    //writeImage(img);
}
