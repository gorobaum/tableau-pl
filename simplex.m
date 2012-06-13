function [ind, x] = simplex(A,b,c,m,n,print)
  [ind, x] = preparesimplex(A,b,c,m,n,print)
endfunction

function [ind, x] = preparesimplex(A,b,c,m,n,print)
  for i = 1:m
    if ( b(i) < 0.0 )
      auxb(i) = -b(i)
      auxA(:,i) = -A(:,i)
    else
      auxb(i) = b(i)
      auxA(:,i) = A(:,i)
    endif
  endfor
  for i = 1:n
    indbase = 0
    auxc(i) = 0.0
  endfor
  for j = 1:m
    indbase = 1
    auxc(j) = 1.0
  endfor

    
  [ind, x, B] = phase1simplex(auxA,auxb,auxc,indbase,m,n,print)
endfunction

function [ind, x, B] = phase1simplex(A,b,c,indbase,m,n,print)
  I = eye(m)
  auxA = [auxA,I]
  auxA = [auxb',A]
  auxc(1) = -(base(c,indbase))'*b
  for i = 2:m+n
    auxc = c(i) - auxc()
  endfor
  auxA = [auxc;auxA]
  
  stop = false
  while (!stop)
    j = 2
    while ( auxA(j,1) > 0.0 )
      j++
    endwhile
    if ( j == length(auxA(:,1)) )
      stop = true
    else
      u = auxA(:,j)
      ratio = inf
      for i = 2:length(u)
        if ( u(i) > 0.0 && u(i)/auxA(i,j) < ratio )
          ratio = u(i)/auxA(i,j)
      endfor
      if ( ratio == inf )
        stop = true
      endif
    endif
  endwhile

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
