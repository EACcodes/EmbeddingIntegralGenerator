/**
 * @file Timer.h
 * @brief Timing utilities
 */

#ifndef TIMER_H
#define	TIMER_H

#include <chrono>
#include <time.h>
#include <string>

/**
 * @class Timer
 * @brief Records wall/cpu times
 */
class Timer {
public:
    
    /**
     * @brief Constructs a new timer
     */
    Timer();
    
    /**
     * @brief Prints out wall/cpu time since the timer was constructed
     * @param message Message to be printed out with the wall/cpu time (usually the name of the code being timed)
     */
    void print_elapsed_time(const std::string& message) const;
    
private:
    /*
     * The C++11 standard gives us a nice high resolution clock for wall times, but not CPU times.
     * For the CPU time we use the classic C clock_t. (Note that this won't work under windows were 
     * clock() doesn't work.)
     */
    typedef std::chrono::high_resolution_clock HighResClock;

    // starting wall and cpu times
    const HighResClock::time_point wall_time_start;
    const clock_t cpu_time_start;
    
};

#endif	/* TIMER_H */

