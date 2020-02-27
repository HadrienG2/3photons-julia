# No include dependency
#
# FIXME: Isn't there a way to spell this out in code???


"Error handling primitives to ease program invariant checking"
module Errors

export check


"""
Simplified variant of standard @assert macro that does not come with a "may be
disabled at various optimization levels" caveat.
"""
macro enforce(expr, msgs...)
    msg_body = isempty(msgs) ? expr : msgs[1]
    msg = string(msg_body)
    return :($(esc(expr)) ? $(nothing) : throw(AssertionError($(msg))))
end

end