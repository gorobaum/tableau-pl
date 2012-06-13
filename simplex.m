function [ind, x] = simplex(A,b,c,m,n,print)
  [ind, x] = preparesimplex(A,b,c,m,n,print)
endfunction

function [ind, x] = preparesimplex(A,b,c,m,n,print)
  for i = 1:m
    if ( b(i) < 0.0 )
      auxb(i) = -b(i)
      auxA(i,:) = -A(i,:)
    else
      auxb(i) = b(i)
      auxA(i,:) = A(i,:)
    endif
  endfor
  for i = 1:n
    auxc(i) = 0.0
  endfor
  for j = (i+1):m+n
    indbase(j-i) = j
    auxc(j) = 1.0
  endfor

  I = eye(m)
  auxA = [auxA,I]
  auxA = [auxb',auxA]
  auxc(1) = -(base(auxc,indbase))'*b
  for i = 2:m+n+1
    auxc(i) = auxc(i-1) - auxc(1)
  endfor
  auxA = [auxc;auxA]
  
  disp("Fase 1")
  [ind, x, B, indb] = runsimplex(auxA,auxb,auxc,indbase,m,m+n,print)
  if (ind == -1 )
    ind = 1
  else
    B = B(:,1:n)
    disp("Fase 2")
    [ind, x, B] = runsimplex(B,b,c,indb,m,n,print)
  endif
endfunction

function [ind, x, B, indb] = runsimplex(A,b,c,indbase,m,n,print)
  
  stop = false
  while (!stop)
    j = 2
    while ( A(j,1) > 0.0 && j < length(A(:,1)))
      j++
    endwhile
    if ( j == length(A(:,1)) )
      stop = true
      ind = 0
    else
      u = A(:,j)
      ratio = inf
      l = 2
      for i = 2:length(u)
        if ( u(i) > 0.0 && u(i)/A(i,j) < ratio )
          ratio = u(i)/A(i,j)
          l = i
        endif
      endfor
      if ( ratio == inf )
        stop = true
        ind = -1
      else
        indbase(l) = j
        for i = 1:length(A(:,j))
          if ( i == l )
            A(i,:) /= A(l,j)
          else
            A(i,:) -= (A(i,j)/A(l,j))*A(l,:)
          endif
        endfor
      endif
    endif
  endwhile
  B = A
  x = A(2:m+1,1)
  indb = indbase
endfunction

function ret = base(vec, b)
  j = 1
  for i = 1:length(vec)
    if ( i == b(j) )
      ret(j) = vec(i)
      j++
    endif
  endfor
endfunction
