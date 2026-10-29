# Go port of the sheaf-core graph subset. `make check` is the gate to run before pushing (remote CI is advisory).
.PHONY: check cross test vet golden clean
check:
	scripts/check.sh
cross:
	CROSS=1 scripts/check.sh
test:
	go test -count=1 ./...
vet:
	go vet ./...
# Regenerate testdata/golden from the Rust reference (needs cargo; Go CI never does this).
golden:
	scripts/gen-golden.sh
clean:
	rm -rf bin
