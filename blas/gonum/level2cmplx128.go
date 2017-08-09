// Copyright ©2017 The gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package gonum

import (
	"gonum.org/v1/gonum/internal/asm/c128"
)

// Zgerc performs the rank-one operation
//  A += alpha * x * y^H
// where A is an m×n dense matrix, alpha is a scalar, x is an m element vector,
// and y is an n element vector.
func (Implementation) Zgerc(m, n int, alpha complex128, x []complex128, incX int, y []complex128, incY int, a []complex128, lda int) {
	// Check inputs
	if m < 0 {
		panic(mLT0)
	}
	if n < 0 {
		panic(nLT0)
	}
	if incX == 0 {
		panic(zeroIncX)
	}
	if incY == 0 {
		panic(zeroIncY)
	}
	if (incX > 0 && (m-1)*incX >= len(x)) || (incX < 0 && (1-m)*incX >= len(x)) {
		panic(badX)
	}
	if (incY > 0 && (n-1)*incY >= len(y)) || (incY < 0 && (1-n)*incY >= len(y)) {
		panic(badY)
	}
	if lda < max(1, n) {
		panic(badLdA)
	}
	if len(a) < (m-1)*lda+n {
		panic("blas: insufficient matrix slice length")
	}

	// Quick return if possible
	if m == 0 || n == 0 || alpha == 0 {
		return
	}

	var kx, jy int
	if incX < 0 {
		kx = (1 - m) * incX
	}
	if incY < 0 {
		jy = (1 - n) * incY
	}
	for j := 0; j < n; j++ {
		if y[jy] != 0 {
			tmp := alpha * conj(y[jy])
			c128.AxpyInc(tmp, x, a[j:], uintptr(n), uintptr(incX), uintptr(lda), uintptr(kx), 0)
		}
		jy += incY
	}
}
