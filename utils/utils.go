package utils

import (
	"fmt"
	"strings"
)

func TrimFloat(x float64, prec int) string {
	s := fmt.Sprintf("%.*f", prec, x)
	s = strings.TrimRight(s, "0")
	s = strings.TrimRight(s, ".")
	return s
}
