#include<Kokkos_Core.hpp>
#include<cstdio>

struct TagA {};
struct TagB {};

struct Foo {
    KOKKOS_INLINE_FUNCTION
    void operator() (const TagA, const Kokkos::TeamPolicy<>::member_type& team) const {
        printf("Greetings from thread %i of team %i with TagA\n",
                team.team_rank(),team.league_rank());
    }
    KOKKOS_INLINE_FUNCTION
    void operator() (const TagB, const Kokkos::TeamPolicy<>::member_type& team) const {
        printf("Greetings from thread %i of team %i with TagB\n",
                team.team_rank(),team.league_rank());
    }
};

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc,argv);

    int N = atoi(argv[1]);

    Foo foo;

    Kokkos::parallel_for(Kokkos::TeamPolicy<TagA>(N,Kokkos::AUTO), foo);
    Kokkos::parallel_for("Loop2", Kokkos::TeamPolicy<TagB>(N,Kokkos::AUTO), foo);

    Kokkos::finalize();
}
