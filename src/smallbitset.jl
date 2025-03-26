mutable struct SmallBitSet <: AbstractSet{Int}
    bits::UInt64

    SmallBitSet() = new(0)
end

"""
    SmallBitSet([itr])

Construct a set represented by a bit string (like Julia's `BitSet`) backed by a
single `UInt64`, in order to achieve a smaller memory size.

!!! warning

    Since this implementation uses a single `UInt64`, storing integers outside
    of the range `1:64` is not supported.
"""
SmallBitSet(itr) = union!(SmallBitSet(), itr)

function Base.empty!(s::SmallBitSet)
    s.bits = 0
    return s
end
Base.empty(::SmallBitSet, ::Type{Int}=Int) = SmallBitSet()
Base.emptymutable(::SmallBitSet, ::Type{Int}=Int) = SmallBitSet()

function Base.copy!(dest::SmallBitSet, src::SmallBitSet)
    dest.bits = src.bits
    return dest
end

Base.copy(s::SmallBitSet) = copy!(SmallBitSet(), s)
Base.copymutable(s::SmallBitSet) = copy(s)

function Base.length(s::SmallBitSet)
    return Base.count_ones(s.bits)
end

function Base.iterate(s::SmallBitSet, state=0)
    index = _next_index(s.bits, state)
    return isnothing(index) ? index : (index, index)
end

function Base.push!(s::SmallBitSet, x::Integer)
    _in_bounds(x) && (s.bits |= _from_index(x))
    return s
end

function Base.push!(s::SmallBitSet, items::Integer...)
    for x in items
        push!(s, x)
    end
    return s
end

function Base.delete!(s::SmallBitSet, x::Integer)
    _in_bounds(x) && (s.bits = s.bits & ~_from_index(x))
    return s
end

function Base.union!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits |= s2.bits
    return s1
end

function Base.union!(s::SmallBitSet, itr)
    for x in itr
        push!(s, x)
    end
    return s
end

function Base.union!(s::SmallBitSet, r::AbstractUnitRange{<:Integer})
    a, b = first(r), last(r)
    a, b = max(a, 1), min(b, 64)
    diff = b - a
    bits = (~0 >>> (63 - diff)) << (a - 1)
    s.bits |= bits
    return s
end

Base.union(s::SmallBitSet, sets...) = union!(copy(s), sets...)

function Base.intersect!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits &= s2.bits
    return s1
end

Base.intersect(s1::SmallBitSet, s2::SmallBitSet) = intersect!(copy(s1), s2)

function Base.setdiff!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits &= ~s2.bits
    return s1
end

function Base.symdiff!(s1::SmallBitSet, s2::SmallBitSet)
    s1.bits = xor(s1.bits, s2.bits)
    return s1
end

Base.issubset(a::SmallBitSet, b::SmallBitSet) = a.bits == a.bits & b.bits
Base.:⊊(a::SmallBitSet, b::SmallBitSet) = a.bits <= b.bits && a.bits != b.bits

@inline Base.in(x::Int, s::SmallBitSet) = _in_bounds(x) ? _get_index(s.bits, x) : false
@inline Base.in(x::Integer, s::SmallBitSet) = _in_bounds(x) ? _get_index(s.bits, x) : false

Base.first(s::SmallBitSet) = _next_index(s.bits, 0)
Base.last(s::SmallBitSet) = _prev_index(s.bits, 63)
Base.minimum(s::SmallBitSet) = first(s)
Base.maximum(s::SmallBitSet) = last(s)
Base.extrema(s::SmallBitSet) = (first(s), last(s))

Base.:(==)(x::SmallBitSet, y::SmallBitSet) = x.bits == y.bits
Base.isequal(x::SmallBitSet, y::SmallBitSet) = isequal(x.bits, y.bits)

@inline _from_index(x::Integer) = unsigned(1 << (x - 1))
@inline _in_bounds(x::Integer) = x > 0 && x <= 64
@inline _get_index(bits::UInt64, idx::Integer) = bits & _from_index(idx) != 0

function _next_index(bits::UInt64, start::Integer)
    start < 0 && (start = 0)
    start > 63 && return nothing
    i = start
    bits >>>= i
    while bits != 0
        if bits & 1 == 1
            return i + 1
        end
        i += 1
        bits >>>= 1
    end
    return nothing
end

function _prev_index(bits::UInt64, start::Integer)
    start < 0 && return nothing
    start > 63 && (start = 63)
    i = start
    bits <<= (63 - i)
    while bits != 0
        if bits & (1 << 63) > 0
            return i + 1
        end
        i -= 1
        bits <<= 1
    end
    return nothing
end

set_to_bits(itr) = reduce(|, (_from_index(x) for x in itr), init=0)

function bits_to_set(bits::Integer)
    s = SmallBitSet()
    s.bits = bits
    return s
end
