program hw2_q2
   use spectral
   
   ! Size
   integer, parameter :: N=8
   integer, parameter :: M=8
   integer, parameter :: L=8

   ! Wavenumbers
   integer, dimension(N) :: k1v
   integer, dimension(M) :: k2v
   integer, dimension(L) :: k3v
   
   ! Function and fourier coefficients
   double complex, dimension(N,M,L) :: f3, f3hat_k 
   double complex, dimension(N,M)   :: f2, f2hat_k 
   double complex, dimension(N)     :: f1, f1hat_k 

   ! Get wavenumbers
   k1v=get_kv(N)
   k2v=get_kv(M)
   k3v=get_kv(L)
   print *, k1v
   ! Initialize coefficients to zero 
   f1hat_k=0.
   f2hat_k=0.
   f3hat_k=0.

   ! Get coefficients that result in real function 
   call initialize_1D(f1hat_k, k1v, N)
   call iFFT_1D(f1, f1hat_k, N)
   ! Output the function
   call write_complex_array(f1hat_k, N, "./outs/q2_1d.txt")
   call write_complex_array(f1, N, "./outs/q2_1d_f.txt")
   
   ! Get coefficients that result in real function 
   call initialize_2D(f2hat_k, k1v, k2v, N, M)
   call iFFT_2D(f2, f2hat_k, N, M)
   ! Output the function
   call write_complex_array_2D(f2hat_k, N, M, "./outs/q2_2d.txt")
   call write_complex_array_2D(f2, N, M, "./outs/q2_2d_f.txt")

   ! Get coefficients that result in real function 
   call initialize_3D(f3hat_k, k1v, k2v, k3v, N, M, L)
   call iFFT_3D(f3, f3hat_k, N, M, L)
   ! Output the function
   call write_complex_array_3D(f3hat_k, N, M, L, "./outs/q2_3d.txt")
   call write_complex_array_3D(f3, N, M, L, "./outs/q2_3d_f.txt")

end program hw2_q2