subroutine getWeightedInt(vals, weights, n, r, res)
implicit none

integer, dimension(n), intent(in) :: vals
double precision, dimension(n), intent(in) :: weights
integer, intent(in) :: n
double precision, intent(in) :: r
integer, intent(out) :: res

integer :: i
double precision :: s

i = 0
s = 0.0

do while (s .le. r)
    do while(abs(weights(i))<= 1e-8)
        i = i+1
    end do
    s = s + weights(i)
    i = i + 1
end do

res = vals(i-1)

end subroutine getWeightedInt

subroutine getWeighted(vals, weights, n, numseq, res)
implicit none

integer, dimension(n), intent(in) :: vals
double precision, dimension(n), intent(in) :: weights
integer, intent(in) :: n
double precision :: r
integer, intent(out) :: res

integer :: i
double precision :: s

r = numseq()

res = getWeightedInt(vals, weights, n, r, res)

end subroutine getWeighted
