function stack_out = normalize_stack(stack)
   energies = zeros(1,size(stack,3))
   stack_out = zeros(size(stack))
   for m = 1 : size(stack,3)
       energies(m) = trapz(trapz(stack(:,:,m))); 
       stack_out(:,:,m) = stack(:,:,m) / energies(m);
   end
