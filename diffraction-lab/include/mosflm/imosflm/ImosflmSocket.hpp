#pragma once

#include "mosflm/control/ControlSession.hpp"

#include <cstdint>
#include <string>

namespace mosflm::imosflm {

/// Minimal iMosflm / TCP session façade (Phase E / E6).
///
/// Full Tcl GUI protocol is not implemented yet; this class owns a
/// :class:`control::ControlSession` and exposes listen/accept/process hooks
/// so a future socket loop can drive the same CONTROL keywords as ``.in`` files.
class ImosflmSocketServer {
public:
    explicit ImosflmSocketServer(std::uint16_t port = 0);

    [[nodiscard]] std::uint16_t port() const { return port_; }

    /// Bind/listen (no-op until real sockets are wired). Returns false if
    /// ``port_ == 0`` (disabled).
    [[nodiscard]] bool listen();

    /// Feed one CONTROL line into the owned session (unit-test / stub path).
    void handle_control_line(const std::string& line);

    [[nodiscard]] control::ControlSession& session() { return session_; }
    [[nodiscard]] const control::ControlSession& session() const {
        return session_;
    }

private:
    std::uint16_t port_{0};
    bool listening_{false};
    control::ControlSession session_;
};

}  // namespace mosflm::imosflm
