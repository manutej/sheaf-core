package sheaf

import (
	"fmt"
	"io"
	"os"
	"strconv"
)

// Usage is the CLI usage line.
const Usage = "usage: sheafgo analyze <spec.json> [--out FILE] [--k K] [--null-reps R] [--seed S] [--top N] [--method shift-invert|polynomial] [--boundary FILE]"

// Main runs the CLI (args excludes the program name) and returns the exit code: 0 ok, 2 input error (one
// "P_INPUT_01 ..." line on stderr), 1 output error. Same flags, codes and messages as the reference `sheaf analyze`.
func Main(args []string, stdout, stderr io.Writer) int {
	fail := func(code int, format string, a ...any) int {
		fmt.Fprintf(stderr, format+"\n", a...)
		return code
	}
	if len(args) < 2 || args[0] != "analyze" {
		return fail(2, "%s", Usage)
	}
	text, err := os.ReadFile(args[1])
	if err != nil {
		return fail(2, "P_INPUT_01 read %s: %v", args[1], err)
	}
	spec, err := Parse(text)
	if err != nil {
		return fail(2, "%v", err)
	}
	var out string
	for i := 2; i < len(args); i += 2 {
		val := ""
		if i+1 < len(args) {
			val = args[i+1]
		}
		badFlag := func(f string) int { return fail(2, "P_INPUT_01 flag %s: bad or missing value %q", f, val) }
		u, uerr := strconv.ParseUint(val, 10, 64)
		switch args[i] {
		case "--out":
			out = val
		case "--k", "--null-reps", "--seed", "--top":
			if uerr != nil {
				return badFlag(args[i])
			}
			switch args[i] {
			case "--k":
				spec.Opts.K = clampInt(u)
			case "--null-reps":
				spec.Opts.NullReps = clampInt(u)
			case "--seed":
				spec.Opts.Seed = u
			case "--top":
				spec.Opts.Top = clampInt(u)
			}
		case "--method":
			if val != "polynomial" && val != "shift-invert" {
				return badFlag("--method")
			}
			spec.Opts.Method = val
		case "--boundary":
			bt, err := os.ReadFile(val)
			if err != nil {
				return fail(2, "P_INPUT_01 read %s: %v", val, err)
			}
			if spec.Boundary, err = ParseBoundary(bt, spec); err != nil {
				return fail(2, "%v", err)
			}
		default:
			return fail(2, "P_INPUT_01 unknown flag %s", args[i])
		}
	}
	rep, err := Analyze(spec)
	if err != nil {
		return fail(2, "%v", err)
	}
	js := rep.JSON()
	if out == "" {
		io.WriteString(stdout, js)
		return 0
	}
	if err := os.WriteFile(out, []byte(js), 0o644); err != nil {
		return fail(1, "P_OUTPUT_01 write %s: %v", out, err)
	}
	fmt.Fprintf(stderr, "wrote %s\n", out)
	return 0
}
