// Command sheafgo is the Go CLI for the graph-sheaf subset:
//
//	sheafgo analyze <spec.json> [--out FILE] [--k K] [--null-reps R] [--seed S] [--top N]
//	                [--method shift-invert|polynomial] [--boundary FILE]
//
// Flags, exit codes and the report format follow the reference `sheaf analyze` (see go/README.md). --method is
// accepted for compatibility; this implementation always uses dense per-component eigendecompositions.
package main

import (
	"os"

	"github.com/manutej/sheaf-core/go/sheaf"
)

func main() { os.Exit(sheaf.Main(os.Args[1:], os.Stdout, os.Stderr)) }
