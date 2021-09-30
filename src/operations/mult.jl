const config = Dict(:multiplication => :fast)

struct MultiplicationType{T} end

get_multiplication_mode() = config

"""
    set_multiplication_mode(multype)

Sets the algorithm used to perform matrix multiplication with interval matrices.

### Input

- `multype` -- symbol describing the algorithm used
   - `:slow` -- uses traditional matrix multiplication algorithm.
   - `:rank1` -- uses rank1 update
   - `:fast` -- computes an enclosure of the matrix product using the midpoint-radius
                notation of the matrix [[RUM10]](@ref).

### Notes

- By default, `:fast` is used.
- Using `fast` is generally significantly faster, but it may return larger intervals,
  especially if midpoint and radius have the same order of magnitude
    (50% overestimate at most) [[RUM99]](@ref).
"""
function set_multiplication_mode(multype)
    multype in (:fast, :slow) || throw(ArgumentError("$multype is not a valid input."))

    type = MultiplicationType{multype}()
    @eval *(A::IntervalMatrix, B::IntervalMatrix) = *($type, A, B)

    @eval *(A::AbstractMatrix, B::IntervalMatrix) = *($type, A, B)

    @eval *(A::IntervalMatrix, B::AbstractMatrix) = *($type, A, B)

    config[:multiplication] = multype
end

*(::MultiplicationType{:slow}, A::IntervalMatrix, B::IntervalMatrix) = IntervalMatrix(A.mat * B.mat)
*(::MultiplicationType{:slow}, A::AbstractMatrix, B::IntervalMatrix) = IntervalMatrix(A * B.mat)
*(::MultiplicationType{:slow}, A::IntervalMatrix, B::AbstractMatrix) = IntervalMatrix(A.mat * B)

function *(::MultiplicationType{:fast},
           A::AbstractIntervalMatrix{Interval{T}},
           B::AbstractIntervalMatrix{Interval{T}}) where {T}

    Ainf = inf(A)
    Asup = sup(A)
    Binf = inf(B)
    Bsup = sup(B)

    mA, mB, R, Csup = setrounding(T, RoundUp) do
        mA = Ainf + 0.5 * (Asup - Ainf)
        mB = Binf + 0.5 * (Bsup - Binf)

        rA = mA - Ainf
        rB = mB - Binf

        R = abs.(mA) * rB + rA * (abs.(mB) + rB)
        Csup = mA * mB + R

        return mA, mB, R, Csup
    end

    Cinf = setrounding(T, RoundDown) do
        mA * mB - R
    end

    return IntervalMatrix(Interval.(Cinf, Csup))
end


function *(::MultiplicationType{:fast},
    A::AbstractMatrix{T},
    B::AbstractIntervalMatrix{Interval{T}}) where {T}

    Binf = inf(B)
    Bsup = sup(B)

    mB, R, Csup = setrounding(T, RoundUp) do
        mB = Binf + 0.5 * (Bsup - Binf)

        rB = mB - Binf

        R = abs.(A) * rB
        Csup = A * mB + R

        return mB, R, Csup
    end

    Cinf = setrounding(T, RoundDown) do
        A * mB - R
    end

    return IntervalMatrix(Interval.(Cinf, Csup))
end

function *(::MultiplicationType{:fast},
    A::AbstractIntervalMatrix{Interval{T}},
    B::AbstractMatrix{T}) where {T}

    Ainf = inf(A)
    Asup = sup(A)

    mA, R, Csup = setrounding(T, RoundUp) do
        mA = Ainf + 0.5 * (Asup - Ainf)

        rA = mA - Ainf

        R = rA * abs.(B)
        Csup = mA * B + R

        return mA, R, Csup
    end

    Cinf = setrounding(T, RoundDown) do
        mA * B - R
    end

    return IntervalMatrix((Interval.(Cinf, Csup)))
end


# function *(::MultiplicationType{:rank1},
#            A::AbstractMatrix{Interval{T}},
#            B::AbstractMatrix{Interval{T}}) where {T<:Real}

#     Ainf = inf.(A)
#     Asup = sup.(A)
#     Binf = inf.(B)
#     Bsup = sup.(B)

#     n = size(A, 2)

#     Csup =  zeros(T, (size(A,1), size(B,2)))
#     Cinf = zeros(T, size(A, 1), size(B, 2))

#     Cinf = setrounding(T, RoundDown) do
#         for i in 1:n
#             Cinf .+= min.(view(Ainf, :, i) * view(Binf, i, :)',
#                           view(Ainf, :, i) * view(Bsup, i, :)',
#                           view(Asup, :, i) * view(Binf, i, :)',
#                           view(Asup, :, i) * view(Bsup, i, :)')
#         end
#         return Cinf
#     end

#     Csup = setrounding(T, RoundUp) do
#         for i in 1:n
#             Csup .+= max.(view(Ainf, :, i) * view(Binf, i, :)',
#                           view(Ainf, :, i) * view(Bsup, i, :)',
#                           view(Asup, :, i) * view(Binf, i, :)',
#                           view(Asup, :, i) * view(Bsup, i, :)')
#         end
#         return Csup
#     end

#     return Interval.(Cinf, Csup)


# end
