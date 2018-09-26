# ARTIFICIAL READS GENERATOR DOCS VERSIÓN PARA LOS DATOS DEL ENCODER

## DESCRIPCIÓN
    Esta versión del ARG genera un archivo de alineamiento .align diferente de la original,
    este archivo tiene una estructura especial para que una máquina de estados extraiga
    facilmente los datos.

## ESTRUCTURA
- B
- C
- '+'   -> Este caracter indica que en la siguiente linea empiezan los reads
- MapPos
- lendesc
- strand
- oper         ->SE REPITEN LENDESC VECES           
- offset       ->SE REPITEN LENDESC VECES   
- baseRef      ->SE REPITEN LENDESC VECES    
- baseRead     ->SE REPITEN LENDESC VECES
- '-' -> Este caracter indica que acabó el read y sigue otro

_NOTA IMPORTANTE: Nunca hacer merge de esta versión al branch master

### DEVELOPERS:
_Juan Camilo Peña Vahos_ @MiloVahos96,

