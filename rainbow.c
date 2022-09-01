#include <math.h>
#include <stdio.h>

void print_rainbow(float shift, float speed) {
  for (int i = 0; i < 100; i++) {
    int r, g, b;
    r = (std::sin(speed * i + 0 + shift) +1) * 127;
            g = (std::sin(speed * i + ((2*M_PI)/3) + shift) +1) * 127;
            b = (std::sin(speed * i + ((4*M_PI)/3)  + shift) +1) * 127;

     printf("r:%4d, g:%4d, b:%4d \n",r,g,b);
   // printf("\033[48;2;%d;%d;%d m  \033[0m", (int)(r), (int)(g), (int)(b));
  }
}
int main() {
  float shift = 0;

  while (1) {
    print_rainbow(shift, 0.1);
  //  printf("\n");
    shift += 0.0001;
  }

  print_rainbow(2, 0.1);

  return 0;
}
