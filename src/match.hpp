#ifndef MATCH_HPP
#define MATCH_HPP

#include "trueskill.hpp"
#include "player.hpp"

namespace trueskill {

class MatchPlayer {
public:
    int id;
    Rating rating;

    MatchPlayer(int playerId, const Rating& initialRating)
        : id(playerId), rating(initialRating) {}

    const Rating& getRating() const {
        return rating;
    }

    void setRating(const Rating& newRating) {
        rating = newRating;
    }
};




} // namespace trueskill

#endif // MATCH_HPP
