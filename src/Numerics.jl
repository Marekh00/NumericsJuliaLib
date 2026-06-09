function RK4(fd,x0,y0,xn)
h=0.01
xpart = x0+h:h:xn
ylist = Array{Float64}(undef,length(xpart))
y = y0
ylist[1]=y
i=1
for n in xpart
k1 = fd(n,y)
k2 = fd(n+h/2,y+(h/2)*k1)
k3 = fd(n+h/2,y+(h/2)*k2)
k4 = fd(n+h,y+h*k3)
y = y + (h/6)*(k1+2k2+2k3+k4)
ylist[i]=y
i+=1
end
return xpart, ylist
end

function IsPrime(n::Integer)::Bool
if n <= 1 return false end
if n <= 3 return true end
if n % 2 == 0 || n % 3 == 0 return false end
range = 5:2:(n/7)
for i in range
if n % i == 0 return false end
end
return true
end

function fac(n)
    if n<0 throw(DomainError(n,"fac is not defined for negative numbers")) end
    if n == 0 || n==1 return 1 end
    return n * fac(n-1)
end

function choose(n,m)::BigInt
    n=big(n);m=big(m)
    return fac(n)/(fac(m)*fac(n-m))
end

function bigpath(tri, i, j , row, col,dp)
    if j == col+1 return 0 end
    if i == row return tri[i,j] end
    
    dp[i,j] = tri[i,j] + max(bigpath(tri,i+1,j,row,col,dp),bigpath(tri,i+1,j+1,row,col,dp))

    return dp[i,j]
end

function Square(M::Matrix)
    return size(M)[1] == size(M)[2]
end

function Transpose(M::Matrix)
    N = zeros(size(M)[2],size(M)[1])
    for i in 1:size(M)[1]
        N[:,i]=M[i,:]
    end
    M = N
    return M
end

function Det(M::Matrix)
    if !Square(M)
        throw(DimensionMismatch("M is not square"))
    end

    sl = size(M)[1]
    if sl == 2
        return M[1,1]*M[2,2]-M[1,2]*M[2,1]
    end

    sum = 0

    for i in 1:sl
        sum += M[1,i]*(Det(M[2:sl,[1:i-1;i+1:sl]])*((-1)^(i+1)))
    end 
    return sum
end

function IdMatrix(n::Int)
    M = zeros(n,n)
    for i in 1:n
        M[i,i]=1
    end
    return M
end

function Invert(M::Matrix)
    if Det(M)==0
        throw(DomainError("M is not invertable"))
    end
    sl = size(M)[1]
    inv = IdMatrix(sl)

    # Mye av problemet med å lage en algoritme som kan invertere en matrise
    # er at diagonalen ikke alltid har kun ikke-null elementer.
    # Vi vil ha ikke-null elementer på diagonalen for å redusere den til identitetsmatrisa
    # ... man kan kanskje bare bearbeide matrisa slik at den har ikke-null elementer FØR man starter ordentlig invertering?

    if sl == 2
        return (1/Det(M))*[M[2,2] -M[1,2] ; -M[2,1] M[1,1]]
    end

    M = Float64.(M)

    for i in 1:sl
        if M[i,i] == 0 # Om diagonal == 0, finn sted hvor det IKKE er det og legg til på original rad
            for j in i+1:sl
                if M[j,i] != 0
                    M[i,:] += M[j,:]
                    inv[i,:] += inv[j,:]
                    break
                end
            end
        end

        val = M[i,i]
        
        M[i,:] *= (1/val)
        inv[i,:] *= (1/val)
        for r in i+1:sl
            coeff = M[r,i]
            M[r,:] -= M[i,:]*coeff
            inv[r,:] -= inv[i,:]*coeff
        end

    end # Etter dette skal alle diagonal-elementene være == 1

    for i in 1:sl-1
        m = sl - i + 1 # Så går vi i revers, velg m-te rad fra bunnen
        for n in 1:m-1 # Trekk fra m-te rad ganger
            coeff = M[n,m]
            M[n,:] -= M[m,:]*coeff
            inv[n,:] -= inv[m,:]*coeff
        end
    end
    return inv

end

function LSM(x::Vector, y::Vector)

    # [x,1]c = y

    L = hcat(x,ones(length(x)))
    K = Transpose(L)
    X = K*L
    X = Invert(X)

    c = X*K*y

    return c

end
