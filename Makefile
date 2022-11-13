CC = gfortran
FLAGS = -o

photo_Halley: ap_photo_Halley.f
	$(CC) $(FLAGS) ap_photo_Halley ap_photo_Halley.f

run:
	./ap_photo_Halley
