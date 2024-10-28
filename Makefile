CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Werror -pedantic -g
TESTFLAGS=-lgtest -lgcov
GCOVFLAGS=--coverage
GCOVDIR = gcov_report
OS = $(shell uname)

all: clean test gcov_report

s21_matrix_oop.a: clean
	$(CXX) $(GCOVFLAGS) -c s21_matrix_oop.cpp
	ar rc s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

test: s21_matrix_oop.a
	$(CXX) tests/*.cpp s21_matrix_oop.a $(CXXFLAGS) $(TESTFLAGS) -o test
	./test

gcov_report: test
	mkdir -p $(GCOVDIR)
	mv *.gcda $(GCOVDIR)
	mv *.gcno $(GCOVDIR)
	lcov --capture --directory $(GCOVDIR) --output-file $(GCOVDIR)/coverage.info
	lcov --remove $(GCOVDIR)/coverage.info '/usr/include/*' --output-file $(GCOVDIR)/coverage_filtered.info
	lcov --remove $(GCOVDIR)/coverage_filtered.info '*/bits/*' --output-file $(GCOVDIR)/coverage_filtered.info
	lcov --remove $(GCOVDIR)/coverage_filtered.info '*/ext/*' --output-file $(GCOVDIR)/coverage_filtered.info
	genhtml $(GCOVDIR)/coverage_filtered.info --output-directory $(GCOVDIR)/html


clang:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -i *.cpp *.h
	rm -rf .clang-format

clang_review:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.cpp *.h
	rm -rf .clang-format

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch $(GCOVDIR) *.html *.css test report *.txt *.dSYM
	find ../../ -type f -name "test_s21_matrix" -delete

check: test
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config s21_matrix_oop.cpp *.h
ifeq ($(OS), Darwin)
	leaks --atExit -- test
else
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test
endif