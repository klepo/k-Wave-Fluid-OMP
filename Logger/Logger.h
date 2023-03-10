/**
 * @file      Logger.h
 *
 * @author    Jiri Jaros \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            jarosjir@fit.vutbr.cz
 *
 * @brief     The header file containing a class responsible for printing out info and error messages
 *            (stdout, and stderr).
 *
 * @version   kspaceFirstOrder 2.17
 *
 * @date      30 April     2017, 11:03 (created) \n
 *            20 February  2019, 14:45 (revised)
 *
 * @copyright Copyright (C) 2019 Jiri Jaros and Bradley Treeby.
 *
 * This file is part of the C++ extension of the [k-Wave Toolbox](http://www.k-wave.org).
 *
 * This file is part of the k-Wave. k-Wave is free software: you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with k-Wave.
 * If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).
 */

#ifndef LOGGER_H
#define LOGGER_H

#include <memory>
#include <iostream>

#include <Logger/OutputMessages.h>
#include <Logger/ErrorMessages.h>

/**
 * @class   Logger
 * @brief   Static class implementing the user interface by info messages.
 * @details StaticClass used for printing out info and error message based on the verbose level. \n
 *          This is a static class.
 */
class Logger
{
  public:
    /**
     * @enum    LogLevel
     * @brief   Log level of the message.
     * @details A enum to specify at which log level the message should be displayed, or the level setd.
     */
    enum class LogLevel
    {
      /// Basic (default) level of verbosity
      kBasic    = 0,
      /// Advanced level of verbosity
      kAdvanced = 1,
      /// Full level of verbosity
      kFull     = 2,
    };

    /// Default constructor is not allowed, static class.
    Logger()              = delete;
    /// Copy constructor is not allowed, static class.
    Logger(const Logger&) = delete;
    /// Destructor is not allowed, static class.
    ~Logger()             = delete;

    /// Operator= is not allowed, static class.
    Logger& operator=(const Logger&) = delete;

    /**
     * @brief Set the log level.
     * @param [in] actualLogLevel - Log level for the logger.
     */
    static void setLevel(const LogLevel actualLogLevel);
    /**
     * Get current log level.
     * @return Current log level.
     */
    static LogLevel getLevel()
    {
      return slogLevel;
    };

    /**
     * @brief   Log desired activity for a given log level, version with string format.
     *
     * @param [in] queryLevel - What level to use.
     * @param [in] format     - Format string.
     * @param [in] args       - Arguments, std::string is not accepted.
     */
    template<typename... Args> static void log(const LogLevel queryLevel, const std::string& format, Args... args)
    {
      if (queryLevel <= Logger::slogLevel)
      {
        std::cout << formatMessage(format, args...);
      }
    } // end of log

    /**
     * @brief Log desired activity for a given log level.
     * @param [in] queryLevel - Log level of the message.
     * @param [in] message    - Message to log.
     */
    static void log(const LogLevel queryLevel, const std::string& message);
    /**
     * @brief Log an error.
     * @param [in] errorMessage - Error message to be printed out.
     */
    static void error(const std::string& errorMessage);
    /**
     * @brief Log an error and terminate the execution.
     * @param [in] errorMessage - error message to be printed to stderr.
     */
    static void errorAndTerminate(const std::string& errorMessage);
    /**
     * @brief Flush output messages.
     * @param [in] queryLevel - Log level of the flush.
     */
    static void flush(const LogLevel queryLevel);

    /**
     * Wrap the line based on delimiters and align it with the rest of the logger output.
     *
     * @param [in] inputString - Input string.
     * @param [in] delimiters  - String of delimiters, every char is a delimiter.
     * @param [in] indentation - Indentation from the beginning.
     * @param [in] lineSize    - Line size.
     * @return Wrapped string.
     *
     * @note The string must not contain tabulator and end-of-line characters.
     */
    static std::string wordWrapString(const std::string& inputString,
      const std::string& delimiters,
      const int indentation = 0,
      const int lineSize    = 65);

    /**
     * @brief   C++-11 replacement for sprintf that works with std::string instead of char*
     * @details The routine was proposed at
     *          http://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
     *          and should work with both Linux and VS 2015.
     *          However it still does not support string in formated arguments
     * @param [in] format - Format string.
     * @param [in] args   - Arguments, std::string is not accepted.
     * @return Formated string.
     */
    template<typename... Args> static std::string formatMessage(const std::string& format, Args... args)
    {
      // when the size is 0, the routine returns the size of the formated string
      size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1; // Extra space for '\0'

      std::unique_ptr<char[]> buf(new char[size]);
      snprintf(buf.get(), size, format.c_str(), args...);
      return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
    }

  private:
    /**
     * @brief Extract a word from a string stream based on delimiters.
     *
     * @param [in,out] textStream - Input text stream.
     * @param [in]     delimiters - List of delimiters as a single string.
     * @return         A word from the string.
     */
    static std::string getWord(std::istringstream& textStream, const std::string& delimiters);

    /// Log level of the logger
    static LogLevel slogLevel;

}; // Logger
//----------------------------------------------------------------------------------------------------------------------

#endif /* LOGGER_H */
