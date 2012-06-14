function [ind, x] = simplex(A,b,c,m,n,print)
  [ind, x] = preparesimplex(A,b,c,m,n,print)
endfunction

function [ind, x] = preparesimplex(A,b,c,m,n,print)
  auxb = [];
  for i = 1:m
    if ( b(i) < 0.0 )
      auxb = [auxb;-b(i)];
      auxA(i,:) = -A(i,:);
    else
      auxb = [auxb;b(i)];
      auxA(i,:) = A(i,:);
    endif
  endfor
  cf1 = [] ;
  for i = 1:n
    cf1 = [cf1;0.0];
  endfor
  for j = (i+1):m+n
    indbase(j-i) = j;
    cf1 = [cf1;1.0];
  endfor
  
  I = eye(m);
  auxA = [auxA,I];
  auxc(1) = -(base(cf1,indbase))'*auxb;
  cf1'-(base(cf1,indbase))'*auxA;
  auxc = [auxc(1),cf1'-(base(cf1,indbase))'*auxA];
  auxA = [auxb,auxA];
  auxA = [auxc;auxA];
  
  disp("Fase 1")
  [ind, x, B, indb] = runsimplex(auxA,auxb,auxc,indbase,m,m+n,print)
  if (ind == -1 )
    ind = 1
  else
    B = B(:,1:n+1)
    disp("Fase 2")
    [ind, x, B] = runsimplex(B,b,c,indb,m,n,print)
  endif
endfunction

function [ind, x, B, indb] = runsimplex(A,b,c,indbase,m,n,print)
  
  stop = false
  while (!stop)
    disp("Novo passo")
    A
    j = 2;
    while ( j <= length(A(1,:)) && A(1,j) >= 0.0 )
      j++;
    endwhile
    if ( j == length(A(1,:))+1 )
      stop = true
      ind = 0
    else
      j
      u = A(:,j)
      ratio = inf;
      for i = 2:length(u)
        if ( u(i) > 0.0 && ratio > A(i,1)/u(i) )
          l = i;
          ratio = A(i,1)/u(i);
        endif
      endfor
      if ( i == length(u)+1 && ratio == inf)
        stop = true
        ind = -1
      else
        l
        indbase(l-1) = j-1
        for i = 1:length(u)
          if ( i != l && u(i) != 0.0 )
            A(i,:) -= (u(i)/u(l))*A(l,:);
          endif
        endfor
        A(l,:) /= u(l);
      endif
    endif
  endwhile
  B = A
  x = A(2:m+1,1)
  indb = indbase
endfunction

function ret = base(vec, b)
  j = 1;
  ret = [];
  for i = 1:length(vec)
    if ( i == b(j) )
      ret = [ret;vec(i)];
      j++;
    endif
  endfor
endfunction
