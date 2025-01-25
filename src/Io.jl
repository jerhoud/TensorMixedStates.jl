export save_state, load_state

"""
    save_state(filename, statename, state)

save the state to disk in a hdf5 file,
several states with different names can be saved in the same file
"""
function save_state(filename::String, statename::String, state::State)
    f = h5open(filename, "cw")
    g = create_group(f, statename)
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

"""
    load_state(filename, statename)

load a state previously saved by `save_sate`
"""
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
