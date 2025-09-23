# GenericOp Simplification

simplify(a::GenericOp) = a

function simplify(a::AnyProd)
    s = simplify.(a.subs)
    simplify_prod(a.coef * prod(prodcoef.(s)), reduce(vcat, prodsubs.(s)))
end

simplify(a::AnySum) =
    simplify_sum(reduce(vcat, sumsubs.(simplify.(a.subs))))

simplify(a::TensorOp{N}) where N = TensorOp{N}(simplify.(a.subs))

simplify(a::PowOp) = simplify_pow(simplfiy(a.arg), a.expo)
simplify(a::ExpOp) = simplify_exp(simplify(a.arg))
simplify(a::DagOp) = simplify_dag(simplify(a.arg))

function simplify(a::Gate)
    sarg = simplify(a.arg)
    simplify_lr(sarg, sarg)
end

function simplify(a::Dissipator)
    sarg = simplify(a.arg)
    darg = simplify_dag(sarg)
    daga = simplify_prod(1, [darg, sarg])
    simplify_lr(sarg, sarg) -
        0.5 * simplify_lr(daga, Identity(daga)) -
        0.5 * simplify_lr(Identity(daga), daga)
end

simplify(a::LeftRight) =
    simplify_lr(simplify(a.left), simplify(a.right))

# IndexepOp simplification

function simplify(a::Evolver)
    sarg = simplify(a.arg)
    simplify_l(sarg) + simplify_r(sarg)
end

simplify(a::Indexed) =
    simplify_ind(simplify(a.op), a.index...)

# GenericOp with indices Simplification

simplify_ind(::Identity, ::Int) = IndexedId
simplify_ind(a::Union{JW_F, Proj, LeftRight{1}}, index) = a(index)
simplify_ind(a::Union{Operator{R, 1}, ExpOp{R, 1}, PowOp{R, 1}}, index) where R =
    a(index)

simplify_ind(a::SumGenOp, index...) =
    simplify_sum(reduce(vcat, (map(x->sumsubs(simplify_ind(x, index...)), a.subs))))

function simplify_ind(a::ProdGenOp, index...)
    s = map(x->simplify_ind(x, index...), a.subs)
    simplify_prod(a.coef * prod(prodcoef.(s)), reduce(vcat, prodsubs.(s)))
end

simplify_ind(a::TensorOp, index...) = simplify_prod(1, tensor_apply(simplify_ind, a, index...))
simplify_ind(::Union{ExpOp, PowOp}, index...) = error("cannot expand multiple site functionals")
simplify_ind(a::Operator, index...) =
    if a.expr isa GenericOp
        a.simplify_ind(a.expr, index...)
    else
        error("cannot simplify an implicit operator")
    end

# helpers

is_involution(::Identity) = true
is_involution(::JW_F) = true
is_involution(a::Operator) = a.type == involution_op
is_involution(::Any) = false

simplify_pow(arg, expo) =
    if expo == 0
        return Identity(arg)
    elseif expo == 1
        return arg
    elseif is_involution(arg)
        simplify_pow(arg, expo - 2 * floor(expo / 2))
    else
        return PowOp(arg, expo)
    end

function simplify_exp(a)
    c = prodcoef(a)
    s = prodsubs(a)
    if length(s) ≠ 1
        return exp(a)
    end
    o = s[1]
    if is_involution(o)
        return cosh(c) * Identity(o) + sinh(c) * o
    end
end 

simplify_dag(a::DagOp) = a.arg
simplify_dag(a::ProdGenOp) = ProdGenOp(conj(a.coef), reverse(simplify_dag.(a.subs)))
simplify_dag(a::SumGenOp) = SumGenOp(simplify_dag.(a.subs))
simplify_dag(a::Operator) =
    if a.type == involution_op || a.type == selfadjoint_op
        a 
    else
        dag(a)
    end

simplify_dag(a) =
    if is_involution(a)
        a 
    else
        dag(a)
    end

simplify_lr()