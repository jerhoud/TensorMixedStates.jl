export save_state, load_state

function save_state(filename::String, statename::String, state::State)
    f = h5open(filename, "cw")
    g = f[statename]
    attrs(g)["show_data"] = state.system.show_data
    if state.type == Pure
        attrs(g)["type"] = 0
        g["pure"] = state.state
        g["mixed"] = MPS(state.system.mixed_sites)
    else
        attrs(g)["type"] = 1
        g["pure"] = MPS(state.system.pure_sites)
        g["mixed"] = state.state
    end
    close(f)
end

function load_state(filename::String, statename::String)
    f = h5open(filename, "r")
    g = f[statename]
    show_data = attrs(g)["show_data"]
    tp = attrs(g)["type"]
    pstate = read(g, "pure", MPS)
    mstate = read(g, "mixed", MPS)
    close(f)
    s = System(siteinds(pstate), siteinds(mstate), show_data)
    if tp == 0
        return State(Pure, s, pstate)
    else
        return State(Mixed, s, mstate)
    end
end
