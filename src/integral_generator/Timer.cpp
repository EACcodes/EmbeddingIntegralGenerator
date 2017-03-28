/* 
 * File:   Timer.cpp
 * Author: dkrisilo
 * 
 * Created on August 13, 2014, 6:10 PM
 */

#include "Timer.h"
#include <iostream>


Timer::Timer(): wall_time_start(HighResClock::now()), cpu_time_start(clock()) {
}

void Timer::print_elapsed_time(const std::string& message) const{
    const HighResClock::time_point wall_time_now = HighResClock::now();
    const clock_t cpu_time_now = clock();

    // Wall time
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(wall_time_now-wall_time_start).count();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(wall_time_now-wall_time_start).count();
    elapsed_milliseconds %= 1000;
    double total_time = elapsed_seconds + elapsed_milliseconds / 1000.0;


    // CPU time
    auto t = cpu_time_now - cpu_time_start;
    auto cpu = float(t)/CLOCKS_PER_SEC;

    std::cout << message << " wall time(s): " << total_time << " cpu time(s): " << cpu << std::endl;

    //std::cout << message << " " << total_time << " seconds" << std::endl;
}
