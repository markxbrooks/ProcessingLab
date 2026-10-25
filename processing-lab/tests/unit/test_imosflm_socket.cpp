#include "mosflm/imosflm/ImosflmSocket.hpp"
#include "test_util.hpp"

using mosflm::imosflm::ImosflmSocketServer;

int main() {
    ImosflmSocketServer srv(0);
    CHECK(!srv.listen());  // disabled when port == 0
    srv.handle_control_line("STRATEGY ANGLE 45");
    CHECK(srv.session().state().strategy_pending);
    // GO runs it; without an orientation and image it says what it needs.
    bool refused = false;
    try {
        (void)srv.session().go();
    } catch (const std::runtime_error& ex) {
        refused = std::string(ex.what()).find("STRATEGY needs") != std::string::npos;
    }
    CHECK(refused);

    ImosflmSocketServer live(9009);
    CHECK(live.listen());
    live.handle_control_line("EXIT");
    CHECK(live.session().exit_requested());

    return mosflm_test::finish("imosflm_socket");
}
