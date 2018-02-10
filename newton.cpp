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
    Image(size_t w, size_t h, std::string n = "")
    : name{n}, width{w}, height{h}
    {
        p.resize(width*height*4);
    }

    void resize(size_t w, size_t h)
    {
        p.resize(w*h*4);
        width = w;
        height = h;
        for (size_t x = 0; x < width; ++x)
        {
            for (size_t y = 0; y < height; ++y)
            {
                setPixel(x, y, (x+y)&1 ? 0 : 0xffffff);
            }
        }
    }

    inline Colour pixel(size_t x, size_t y) const
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
    double left;
    double top;
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
    double xStep() const
    {
        return width / static_cast<double>(pixelWidth);
    }
    double yStep() const
    {
        return height / static_cast<double>(pixelHeight);
    }
};

void calculateNewtonChunk(int xBegin, int xEnd, int yBegin, int yEnd,
                          FractalInfo info, Image& img,
                          const std::vector<double>& aaPos,
                          int* done, const bool& shouldExit)
{
    double xs = info.xStep();
    double ys = info.yStep();
    double xb = info.left;
    double yb = info.top;
    size_t it = 0;
    for (int x = xBegin; x < xEnd && !shouldExit; ++x)
    {
        for (int y = yBegin; y < yEnd; ++y)
        {
            ++it;
            double cx = xs * x + xb;
            double cy = ys * y + yb;
            img.setPixel(x, y, newtonAA(xs, ys, info.samples, aaPos, info.chaos,
                                        cx, cy, info.exponent,
                                        info.maxIterations, 1e-16));

        }
    }
    *done = shouldExit ? 2 : 1;
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
    info.height = 4.;
    info.pixelWidth = 800;
    info.pixelHeight = 600;
    info.exponent = 10;
    info.samples = 51;
    info.maxIterations = 96;
    info.chaos.real = -8.;
    info.chaos.imag = 0.;
    info.recalculateWidth();
    info.left = -0.5 * info.width;
    info.top = -0.5 * info.height;
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

    void requestStop()
    {
        shouldExit = true;
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

        shouldExit = false;
        if (info.pixelWidth != img.width || info.pixelHeight != img.height)
        {
            std::cout << "Resizing from " << img.width << "x" << img.height << " to " << info.pixelWidth << "x" << info.pixelHeight << std::endl;
            img.resize(info.pixelWidth, info.pixelHeight);
        }

        int yBegin = 0;
        auto maxThreads = std::thread::hardware_concurrency();
        done.resize(maxThreads+1, false);

        int chunkHeight = info.pixelHeight / maxThreads;
        auto positions = getAAPos(info.samples);
        for (unsigned i = 0; i < maxThreads; ++i)
        {
            threads.emplace_back(calculateNewtonChunk, 0, info.pixelWidth,
                                 yBegin, yBegin + chunkHeight, info, std::ref(img),
                                 positions, &done[i], std::ref(shouldExit));
            yBegin += chunkHeight;
        }

        threads.emplace_back(calculateNewtonChunk, 0, info.pixelWidth,
                             yBegin, info.pixelHeight, info, std::ref(img),
                             positions, &done[done.size()-1],
                             std::ref(shouldExit));
    }

private:
    std::vector<int> done;
    std::vector<std::thread> threads;
    bool shouldExit = false;
};

Image renderPreliminary(FractalInfo info, int rx = 2, int ry = 2)
{
    info.samples = 1;
    info.width  *=  1 << rx;
    info.height *=  1 << ry;
    auto scaledWidth = (info.pixelWidth >> rx) + 1;
    auto scaledHeight = (info.pixelHeight >> ry) + 1;
    Image small(scaledWidth, scaledHeight);
    int dummyDone = 0;
    for (size_t x = 0; x < scaledWidth; ++x)
    {
        for (size_t y = 0; y < scaledHeight; ++y)
        {
            small.setPixel(x, y, x+y);
        }
    }
    auto pos = getAAPos(info.samples);
    calculateNewtonChunk(0, scaledWidth, 0, scaledHeight,
                         info, small, pos,
                         &dummyDone, false);

    Image upscaled(info.pixelWidth, info.pixelHeight);
    for (size_t x = 0; x < info.pixelWidth; ++x)
    {
        for (size_t y = 0; y < info.pixelHeight; ++y)
        {
            upscaled.setPixel(x, y, small.pixel(x>>rx, y>>ry));
        }
    }
    return upscaled;
}

void zoom(FractalInfo& info, Complex<double> center, double amount)
{
    amount = 1. / amount;
    auto newWidth = amount * info.width;
    auto newHeight = amount * info.height;
    info.left += center.real * info.width * (1. - amount);
    info.top += center.imag * info.height * (1. - amount);
    info.width = newWidth;
    info.height = newHeight;
}

void getDragRectangle(Complex<double> dragBegin, Complex<double> dragEnd,
                      double& x0, double& x1, double& y0, double& y1)
{
    x0 = dragBegin.real;
    x1 = dragEnd.real;
    y0 = dragBegin.imag;
    y1 = dragEnd.imag;

    double dxMin = std::min(dragBegin.real, dragEnd.real);
    double dxMax = std::max(dragBegin.real, dragEnd.real);
    double dyMin = std::min(dragBegin.imag, dragEnd.imag);
    double dyMax = std::max(dragBegin.imag, dragEnd.imag);

    double xLen = dxMax - dxMin;
    double yLen = dyMax - dyMin;

    bool xReverse = x0 > x1;
    bool yReverse = y0 > y1;

    if (xLen > yLen)
    {
        if (xReverse) std::swap(x0, x1);
        if (yReverse) y1 = y0 - (x1 - x0);
        else          y1 = y0 + (x1 - x0);
    }

    else
    {
        if (yReverse) std::swap(y0, y1);
        if (xReverse) x1 = x0 - (y1 - y0);
        else          x1 = x0 + (y1 - y0);
    }
}

std::string trimLeading(std::string s, char c)
{
    size_t idx = 0;
    while (idx < s.size() && s[idx] == c) ++idx;
    if (idx == s.size()) return "";
    return s.substr(idx);
}

template <typename T>
void readFromArgs(T& val, int& c, int argc, char* argv[])
{
    ++c;
    if (c >= argc)
    {
        // Yes, I could throw an exception and let calling code deal with it but
        // this is simpler and this function is only called once in order to
        // read program arguments.
        std::cerr << "Invalid argument.\n";
        std::exit(1);
    }
    std::stringstream ss(argv[c]);
    ss >> val;
    if (ss.fail() || ss.bad())
    {
        std::cerr << "Unable to parse '" << argv[c] << "'.\n";
        std::exit(2);
    }
}

FractalInfo extractInfo(int argc, char* argv[])
{
    FractalInfo info = defaultInfo();
    bool getFromStdin = argc == 1;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        arg = trimLeading(arg, '-');
        if (arg == "stdin") getFromStdin = true;
        else if (arg == "left") readFromArgs(info.left, i, argc, argv);
        else if (arg == "top") readFromArgs(info.top, i, argc, argv);
        else if (arg == "width") readFromArgs(info.width, i, argc, argv);
        else if (arg == "height") readFromArgs(info.height, i, argc, argv);
        else if (arg == "exponent") readFromArgs(info.exponent, i, argc, argv);
        else if (arg == "samples") readFromArgs(info.samples, i, argc, argv);
        else if (arg == "iterations" || arg == "maxIterations")
        {
            readFromArgs(info.maxIterations, i, argc, argv);
        }
        else if (arg == "chaos")
        {
            readFromArgs(info.chaos.real, i, argc, argv);
            readFromArgs(info.chaos.imag, i, argc, argv);
        }
        else if (arg == "r" || arg == "real")
        {
            readFromArgs(info.chaos.real, i, argc, argv);
        }
        else if (arg == "i" || arg == "imag")
        {
            readFromArgs(info.chaos.imag, i, argc, argv);
        }
        else if (arg == "default"); // Do nothing; default arguments are used.
        else
        {
            std::cerr << "Unrecognised argument: " << argv[i] << "\n";
            std::exit(3);
        }
    }
    if (getFromStdin) info = getInfoFromStdin();
    return info;
}

int main(int argc, char* argv[])
{
    std::ios_base::sync_with_stdio(false);
    ImageRenderer ir;
    ir.info = extractInfo(argc, argv);

    sf::RenderWindow wnd(sf::VideoMode(ir.info.pixelWidth, ir.info.pixelHeight),
                         "Newton fractal generator");
    wnd.setVerticalSyncEnabled(true);

    auto lastSize = wnd.getSize();

    bool needsRender = true;

    bool isDragging = false;
    Complex<double> dragBegin{0., 0.}, dragEnd{0., 0.};

    while (wnd.isOpen())
    {
        sf::Event evt;
        while (wnd.pollEvent(evt))
        {
            if (evt.type == sf::Event::Closed)
            {
                wnd.close();
            }
            if (evt.type == sf::Event::KeyPressed)
            {
                needsRender = true;
                switch (evt.key.code)
                {
                case sf::Keyboard::Left:
                    ir.info.left -= 0.025 * ir.info.width;
                    break;
                case sf::Keyboard::Right:
                    ir.info.left += 0.025 * ir.info.width;
                    break;
                case sf::Keyboard::Up:
                    ir.info.top -= 0.025 * ir.info.height;
                    break;
                case sf::Keyboard::Down:
                    ir.info.top += 0.025 * ir.info.height;
                    break;
                case sf::Keyboard::I:
                    zoom(ir.info, {0.5, 0.5}, 1.25);
                    break;
                case sf::Keyboard::O:
                    zoom(ir.info, {0.5, 0.5}, 0.8);
                    break;
                case sf::Keyboard::Escape:
                    wnd.close();
                    break;
                default:
                    needsRender = false;
                }
            }
            if (evt.type == sf::Event::MouseButtonPressed)
            {
                double mx = evt.mouseButton.x / (double)wnd.getSize().x;
                double my = evt.mouseButton.y / (double)wnd.getSize().y;
                if (evt.mouseButton.button == sf::Mouse::Left)
                {
                    isDragging = true;
                    dragBegin = {mx, my};
                    dragEnd = dragBegin;
                }
                if (evt.mouseButton.button == sf::Mouse::Right)
                {
                    zoom(ir.info, {mx, my}, 0.8);
                    needsRender = true;
                }
            }
            if (evt.type == sf::Event::MouseMoved
                && sf::Mouse::isButtonPressed(sf::Mouse::Left))
            {
                double mx = evt.mouseMove.x / (double)wnd.getSize().x;
                double my = evt.mouseMove.y / (double)wnd.getSize().y;
                dragEnd = {mx, my};
            }
            if (evt.type == sf::Event::MouseButtonReleased
                && evt.mouseButton.button == sf::Mouse::Left)
            {
                isDragging = false;
                double mx = evt.mouseButton.x / (double)wnd.getSize().x;
                double my = evt.mouseButton.y / (double)wnd.getSize().y;
                dragEnd = {mx, my};
                double x0, x1, y0, y1;
                getDragRectangle(dragBegin, dragEnd,
                                 x0, x1, y0, y1);
                auto newWidth = x1-x0;
                auto newHeight = y1-y0;
                // Bad rectangle, we just continue processing events instead of
                // zooming.
                if (newWidth < 1e-6 || newHeight < 1e-6) continue;
                ir.info.left += x0 * ir.info.width;
                ir.info.top += y0 * ir.info.height;
                ir.info.width *= newWidth;
                ir.info.height *= newHeight;
                needsRender = true;
            }
        }
        wnd.clear(sf::Color::Black);
        auto sz = wnd.getSize();
        wnd.setView(sf::View(sf::FloatRect(0, 0, sz.x, sz.y)));

        if (lastSize != sz)
        {
            lastSize = sz;
            needsRender = true;
        }


        if (needsRender)
        {
            FractalInfo thumbInfo = ir.info;
            thumbInfo.pixelWidth = sz.x;
            thumbInfo.pixelHeight = sz.y;
            thumbInfo.recalculateWidth();
            auto thumb = renderPreliminary(thumbInfo, 3, 3);
            if (ir.isReady())
            {
                lastSize = sz;
                ir.info.pixelWidth = sz.x;
                ir.info.pixelHeight = sz.y;
                ir.info.recalculateWidth();

                ir.img = thumb;

                ir.calculate();
                needsRender = false;
            }
            else
            {
                ir.requestStop();
            }
        }
        sf::Texture texture;
        texture.create(ir.img.width, ir.img.height);
        texture.update(ir.img.p.data());
        sf::Sprite sprite;
        sprite.setTexture(texture);

        wnd.draw(sprite);
        if (isDragging)
        {
            sf::VertexArray rect(sf::LinesStrip, 5);
            double x0, x1, y0, y1;
            double aspect = (double)sz.y / sz.x;
            getDragRectangle(dragBegin, dragEnd,
                             x0, x1, y0, y1);

            rect[0].position = sf::Vector2f(x0, y0);
            rect[1].position = sf::Vector2f(x0, y1);
            rect[3].position = sf::Vector2f(x1, y0);
            rect[2].position = sf::Vector2f(x1, y1);
            for (size_t i = 0; i < 4; ++i)
            {
                rect[i].color = sf::Color::Red;
                rect[i].position.x *= sz.x;
                rect[i].position.y *= sz.y;
            }
            rect[4] = rect[0];
            wnd.draw(rect);
        }
        wnd.display();
    }
    ir.requestStop();
}
