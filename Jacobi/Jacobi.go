package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

type Matrix_CSR struct {
	n, nnz    int64
	IA, JA    []int64
	AIJ, Diag []float64
}

func reading() *Matrix_CSR {
	file, err := os.Open("/mnt/c/Users/dns/Documents/GitHub/Jacobi/matrix/matrix1_n5.mycsr")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	Matrix := new(Matrix_CSR)
	IA := false
	JA := false
	AIJ := false

	rd := bufio.NewReader(file)
	for {
		line, err := rd.ReadString('\n')

		if err != nil {
			if err == io.EOF {
				break
			}
			log.Fatalf("read file line error: %v", err)
			return nil
		}

		splitted_line := strings.Split(line, " ")
		for _, value := range splitted_line {
			if splitted_line[0] == "n" {
				n, err := strconv.ParseUint(splitted_line[2], 10, 64)
				if err != nil {
					fmt.Println(err)
				}
				Matrix.n = int64(n)
				nnz, err := strconv.ParseUint(splitted_line[5], 10, 64)
				if err != nil {
					fmt.Println(err)
				}
				Matrix.nnz = int64(nnz)
				break
			} else if splitted_line[0] == "VECTOR" {
				if splitted_line[1] == "IA" {
					IA = true
					JA = false
					AIJ = false
				}
			} else {
				if IA == true {
					n, err := strconv.ParseUint(value, 10, 64)
					if err != nil {
						fmt.Println(err)
					}
					Matrix.IA = append(Matrix.IA, int64(n))
				} else if JA == true {
					n, err := strconv.ParseUint(value, 10, 64)
					if err != nil {
						fmt.Println(err)
					}
					Matrix.JA = append(Matrix.JA, int64(n))
				} else if AIJ == true {
					n, err := strconv.ParseUint(value, 10, 64)
					if err != nil {
						fmt.Println(err)
					}
					Matrix.AIJ = append(Matrix.AIJ, float64(n))
				}
			}
		}
	}

	return Matrix
}

func main() {
	Matrix := reading()
	fmt.Println(&Matrix)
}
