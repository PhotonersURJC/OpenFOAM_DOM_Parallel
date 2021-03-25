Cambiada quadrature

Separadas calculate de calculateBand para evitar ifs
Actualización de Qr va a updateBandG y updateG, para calcular entre iteraciones sin necesidad
Pasan a llamarse updateQrAndG y updateBandQrAndG

Creo que las bandas no funcionan porque los campos por bandas de la cuadratura tienen tamaño 1.

No tiene paralelización