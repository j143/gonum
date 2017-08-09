package testblas

import "testing"

type Zgercer interface {
	Zgerc(m, n int, alpha complex128, x []complex128, incX int, y []complex128, incY int, a []complex128, lda int)
}

func ZgercTest(t *testing.T, impl Zgercer) {
	for tc, test := range []struct {
		m, n  int
		alpha complex128
		x, y  []complex128
		a     []complex128

		want []complex128
	}{} {
	}
}
