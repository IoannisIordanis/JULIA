using LinearAlgebra
function compiv(A)
  #ολική οδήγηση Gauss
  #Γιάννης Ιορδάνης
  A=Float64.(A)
  n=size(A,1)
  p=Array{Int64}(ones(n));for i=1:n;p[i]=i;end
  q=Array{Int64}(ones(n));for i=1:n;q[i]=i;end
  for k= 1:n-1
    #-------------------------------------------------
    #εντοπισμός τιμής μέγιστου (απόλυτα) σε υποπίνακα
    local ma=maximum(abs.(A[k:n,k:n]))
    #------------------------------------------------
    #εντοπισμός θέσης μέγιστου (απόλυτα) σε υποπίνακα
    local pos = findfirst( x -> x == ma, abs.(A[k:n,k:n]))
    local i=pos[1]; local j=pos[2]
    local i = i+k-1; local j = j+k-1
    #------------------------------------------------
    #αντιμετάθεση γραμμών  και στηλών
    A[[k i], :] = A[[i k], :]
    A[:,[k j]] = A[:,[j k]]
    #--------------------------------------
    # παραγωγή μεταθετικών πινάκων P και Q
    p[[k i]] = p[[i k]]
    q[[k j]] = q[[j k]]
    #----------------------------------------
    #υπολογισμός πολλαπλασιαστών και απαλοιφή
    if A[k, k] != 0.0
      A[k+1:n, k] = A[k+1:n, k]/A[k, k]
      A[k+1:n, k+1:n] = A[k+1:n, k+1:n] - A[k+1:n, k]*(A[k, k+1:n])';
    end
    #----------------------------------------------
  end
  #----------------------------------------------
  U=UpperTriangular(A)
  L=UnitLowerTriangular(A)
  P1=Diagonal(ones(n,n))+zeros(n,n); P=P1[p,:]
  Q1=Diagonal(ones(n,n))+zeros(n,n); Q=Q1[:,q]

  return P,Q,L,U
end