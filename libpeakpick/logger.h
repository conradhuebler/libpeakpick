/*
 * <Simple logging system>
 * Copyright (C) 2024  Conrad Hübler <Conrad.Huebler@gmx.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <iostream>
#include <sstream>
#include <string>

namespace PeakPick {

enum class LogLevel {
    Debug = 0,
    Info = 1,
    Warning = 2,
    Error = 3,
    None = 4
};

class Logger {
private:
    static LogLevel current_level;
    static bool use_colors;

public:
    static void setLevel(LogLevel level)
    {
        current_level = level;
    }

    static LogLevel getLevel()
    {
        return current_level;
    }

    static void setColors(bool enable)
    {
        use_colors = enable;
    }

    static void log(LogLevel level, const std::string& message, const std::string& context = "")
    {
        if (level < current_level)
            return;

        std::string prefix;
        std::string color_code;
        std::string reset_code = use_colors ? "\033[0m" : "";

        switch (level) {
        case LogLevel::Debug:
            prefix = "[DEBUG]";
            color_code = use_colors ? "\033[36m" : ""; // Cyan
            break;
        case LogLevel::Info:
            prefix = "[INFO]";
            color_code = use_colors ? "\033[32m" : ""; // Green
            break;
        case LogLevel::Warning:
            prefix = "[WARNING]";
            color_code = use_colors ? "\033[33m" : ""; // Yellow
            break;
        case LogLevel::Error:
            prefix = "[ERROR]";
            color_code = use_colors ? "\033[31m" : ""; // Red
            break;
        default:
            return;
        }

        std::ostream& out = (level >= LogLevel::Warning) ? std::cerr : std::cout;

        out << color_code << prefix;
        if (!context.empty()) {
            out << " [" << context << "]";
        }
        out << " " << message << reset_code << std::endl;
    }

    // Convenience methods
    static void debug(const std::string& message, const std::string& context = "")
    {
        log(LogLevel::Debug, message, context);
    }

    static void info(const std::string& message, const std::string& context = "")
    {
        log(LogLevel::Info, message, context);
    }

    static void warning(const std::string& message, const std::string& context = "")
    {
        log(LogLevel::Warning, message, context);
    }

    static void error(const std::string& message, const std::string& context = "")
    {
        log(LogLevel::Error, message, context);
    }
};

// Initialize static members
LogLevel Logger::current_level = LogLevel::Info;
bool Logger::use_colors = true;

// Macros for convenient logging with context
#define LOG_DEBUG(msg) PeakPick::Logger::debug(msg, __FUNCTION__)
#define LOG_INFO(msg) PeakPick::Logger::info(msg, __FUNCTION__)
#define LOG_WARNING(msg) PeakPick::Logger::warning(msg, __FUNCTION__)
#define LOG_ERROR(msg) PeakPick::Logger::error(msg, __FUNCTION__)

} // namespace PeakPick
