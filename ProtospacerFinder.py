
#Módulos######################################################################
##############################################################################
import os
import sys
import datetime
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

#Funciones####################################################################
##############################################################################
def linear_FASTA (fichero):
    """
    (str)->{dic}
    Esta funcion toma como entrada la direccion de un fichero con formato FASTA,
    lo abre y devuelve un diccionario en el que las claves son las cabeceras y
    los valores las secuencias
    """
    dic = {}
    head = ""
    for line in open(fichero):
        line = line.rstrip()
        if ">" in line:
            head = line[1:]
            seq = ""
        else:
            line = line.upper()
            seq = seq+line
            dic[head] = seq
    return dic

def complementaria (seq):
    """
    (str)->(str)
    Esta funcion toma una secuencia de ADN, una cadena de texto, y devuelve la
    complementaria de la inversa de la misma
    """
    bases_complementarias = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    newseq = ''
    for base in seq:
        newseq = bases_complementarias[base] + newseq
    return newseq

def do_BLAST (db, sp):
    """
    (str),(str)
    Esta función toma como entrada dos cadenas de texto correspondientes al fichero
    FASTA con los datos genómicos/metagenómicos y la secuencia para ejecutar un BLASTn.
    Antes de hacer el BLAST, crea la base de datos necesaria con el fichero FASTA
    """
    db_cmd = "makeblastdb -in '"+db+"' -dbtype nucl"
    os.system(db_cmd)   
    blast_cmd = "blastn -query "+sp+" -db "+db+" -gapopen 10 -gapextend 2 -reward 1 -penalty -1 -evalue 1 -word_size 7 -out temp.xml -outfmt 5"
    os.system(blast_cmd)

def extract_seq (name):
    """
    (str)->(str)
    Funcion que elimina el ">" del identificador de una secuencia de un fichero formato
    FASTA
    """
    nlist = name.split()
    del nlist[0]
    seq_name = " ".join(nlist)
    return seq_name   

def pam (seq, pos_i,pos_f):
    """
    (str),(int),(int)->(str),(str)
    Funcion que toma como entrada una secuencia perteneciente a un genoma y las
    posiciones de inicio y final del protoespaciador, y devuelve las secuencias
    PAM 5' (aguas arriba) y 3` (aguas abajo) de dicho protoespaciador
    """
    if pos_i < pos_f: #cuando el protoespaciador esta en sentido 5'-3'
        pam5 = seq[pos_i-4:pos_i-1]
        pam3 = seq[pos_f:pos_f+3]
    else: #cuando el protoespaciador esta en sentido 3'-5'
        pam5 = complementaria(seq[pos_i:pos_i+3])
        pam3 = complementaria(seq[pos_f-4:pos_f-1])
    return pam5, pam3

#Inputs#######################################################################
##############################################################################
db_dir = input("Directorio con los genomas/bases de datos que analizar: ")
db_dir = os.path.abspath(db_dir)
sp = input("Fichero con los espaciadores: ")
sp_dire = os.path.abspath(sp)
d = int(input("Número máximo de discrepancias en la diana: "))
outname = input("Nombre del fichero de salida: ")
save = (input("¿Dese guardar las secuencias?[y/n]: ")).lower()

#Ejecucion####################################################################
##############################################################################
start_time = datetime.datetime.now()

#Se crean y se abren los ficheros de salida del programa
file_name1 = outname + ".txt"
file_name2 = outname + "_sequences.txt"
f1 = open(file_name1, "w+")
if save == "y":
    f2 = open(file_name2, "w+")

#Se recorre cada fichero del directorio indicado
for db in os.listdir(db_dir):
    if db.endswith(".fasta"):
        db_dire = db_dir + "/" + db
        
         #Se crea un diccionario cabecera:secuencia
        genome_dic = linear_FASTA(db_dire)
        #Se crea una base de datos db_blast y se realiza el BLAST
        do_BLAST(db_dire, sp_dire)
        
        print("Blast done")
        
        #La informacion resultante del BLAST, formato xml, se convierte en un objeto blast_records
        result_handle = open("temp.xml")
        blast_records = NCBIXML.parse(result_handle)
        
        #Se extrae la información del objeto blast_records y se imprime en el fichero de salida
        for blast_record in blast_records:
            sp_len = blast_record.query_letters
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    hits = hsp.identities
                    if (sp_len - hits <= d):
                        
                        #Printing alignment
                        
                        print("Genome file:", db, file=f1)
                        print("Query spacer:", blast_record.query, "(Spacer length:",sp_len, "pb)", file=f1)
                        print("Target sequence:", alignment.title, "\n", file=f1)
                        sp_s = hsp.query_start
                        sp_e = hsp.query_end
                        seq_s = hsp.sbjct_start
                        seq_e= hsp.sbjct_end
                        if len(str(seq_s))>=7:
                            print(sp_s,"\t\t", hsp.query, "\t\t", sp_e, file=f1)
                            print("\t\t",hsp.match, file=f1)
                            print(seq_s,"\t", hsp.sbjct, "\t\t", seq_e, "\n", file=f1)
                        else:
                            print(sp_s,"\t\t", hsp.query, "\t\t", sp_e, file=f1)
                            print("\t\t",hsp.match, file=f1)
                            print(seq_s,"\t\t", hsp.sbjct, "\t\t", seq_e, "\n", file=f1)
                            
                        #Printing PAMs: se obtienen las secuencias correespondientes a las cabeceras y se obtienen las secuencias PAM
                            
                        if sp_s != 1:
                            seq_i = seq_s - (sp_s - 1)
                        else:
                            seq_i = seq_s
                        if sp_e != sp_len:
                            seq_f = seq_e + (sp_len - sp_e)
                        else:
                            seq_f = seq_e
                        seq = extract_seq(alignment.title)
                        pam5,pam3 = pam(genome_dic[seq],seq_i,seq_f)
                        print("PAM 5': ", pam5, file=f1)
                        print("PAM 3': ", pam3, file=f1)
                        print("*************************************************************************\n", file=f1)
                        
                        #Printing the subject sequence
                        
                        if save == "y":
                            print(">",seq, file=f2)
                            print(genome_dic[seq], file=f2)
                        
        os.remove("temp.xml") 
        
end_time = datetime.datetime.now()
print(end_time-start_time)
