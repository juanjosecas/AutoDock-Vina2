#ifndef VINA_CONSOLE_H
#define VINA_CONSOLE_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#ifdef _WIN32
#include <io.h>
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace vina_console {
inline bool interactive() {
    const char* term = std::getenv("TERM");
    if (term && std::strcmp(term, "dumb") == 0) return false;
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

inline bool supports_color() {
    if (std::getenv("NO_COLOR") || !interactive()) return false;
#ifdef _WIN32
    const HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    if (!GetConsoleMode(handle, &mode)) return false;
    // ENABLE_VIRTUAL_TERMINAL_PROCESSING, including older Windows SDKs.
    return SetConsoleMode(handle, mode | 0x0004) != 0;
#else
    return true;
#endif
}

inline std::string color(const std::string& text, const char* code) {
    static const bool enabled = supports_color();
    return enabled ? std::string(code) + text + "\033[0m" : text;
}
inline std::string heading(const std::string& text) { return color(text, "\033[1;36m"); }
inline std::string success(const std::string& text) { return color(text, "\033[1;32m"); }
}
#endif
