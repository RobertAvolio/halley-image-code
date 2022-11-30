CC = gfortran
FLAGS = -o

photo_Halley: ap_photo_Halley.f
	$(CC) $(FLAGS) ap_photo_Halley ap_photo_Halley.f

auto_Halley: ap_photo_Halley_auto.f
	$(CC) $(FLAGS) ap_photo_Halley_auto ap_photo_Halley_auto.f

run:
	./ap_photo_Halley_auto
