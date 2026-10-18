#include "mosflm/imosflm/ImosflmSocket.hpp"

#include <stdexcept>

namespace mosflm::imosflm {

ImosflmSocketServer::ImosflmSocketServer(std::uint16_t port) : port_(port) {}

bool ImosflmSocketServer::listen() {
    if (port_ == 0) {
        return false;
    }
    // Real BSD sockets / iMosflm binary protocol TBD (Phase E / E6).
    // Binding is deferred so builds stay portable without linking network code
    // until the protocol dump is captured from Fortran open_socket.
    listening_ = true;
    return true;
}

void ImosflmSocketServer::handle_control_line(const std::string& line) {
    if (!listening_ && port_ != 0) {
        throw std::runtime_error("ImosflmSocketServer: listen() not called");
    }
    session_.parse_line(line);
}

}  // namespace mosflm::imosflm
