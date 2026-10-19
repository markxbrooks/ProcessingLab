#include "mosflm/imosflm/ImosflmSocket.hpp"
#include "test_util.hpp"

using mosflm::imosflm::ImosflmSocketServer;

int main() {
    ImosflmSocketServer srv(0);
    CHECK(!srv.listen());  // disabled when port == 0
    srv.handle_control_line("STRATEGY ANGLE 45");
    CHECK(srv.session().state().strategy_pending);
    const auto r = srv.session().go();
    CHECK(r.ran_strategy);

    ImosflmSocketServer live(9009);
    CHECK(live.listen());
    live.handle_control_line("EXIT");
    CHECK(live.session().exit_requested());

    return mosflm_test::finish("imosflm_socket");
}
