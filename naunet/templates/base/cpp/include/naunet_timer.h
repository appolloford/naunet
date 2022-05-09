#ifndef __NAUNET_TIMER_H__
#define __NAUNET_TIMER_H__

#include <chrono>

class Timer {
   private:
    bool m_status;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
    std::chrono::duration<double> m_duration;

   public:
    void start() {
        m_start  = std::chrono::high_resolution_clock::now();
        m_status = true;
    };
    void stop() {
        m_end    = std::chrono::high_resolution_clock::now();
        m_status = false;
    };
    void restart() {
        m_start  = std::chrono::high_resolution_clock::now();
        m_status = false;
    }
    double elapsed() {
        // return
        // std::chrono::duration_cast<std::chrono::microseconds>(m_end-m_start).count();
        m_duration = m_end - m_start;
        return m_duration.count();
    };
};

#endif
