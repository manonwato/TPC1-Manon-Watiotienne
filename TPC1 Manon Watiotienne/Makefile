CC = gcc
# CC = mingw32-gcc
# CFLAGS = -Wall -Wextra
# CFLAGS = ""

SRC = blas1_static.c blas2_static.c jacobi_dynamic.c jacobi_static.c utilities.c main.c GaussSeidel_static.c GaussSeidel_dynamic.c SOR_static.c SOR_dynamic.c

TARGET = tpc1

$(TARGET): $(SRC)
	$(CC) $(SRC) -o $@
# 	$(CC) $(CFLAGS) $(SRC) -o $@

clean: clean
	rm -f $(TARGET)

# re: fclean all