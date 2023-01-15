#imports
import os
from urllib import request
from urllib import error
from Bio import SeqIO

# import DLA_Ranker as RNK

pdbMatrix = []
fastaMatrix = []

# Defining the function for downloading the Fasta files corresponding to the pdb files
def download_url(url, root, filename):
    """Download a file from a url and place it in root.
    Args:
        url (str): URL to download file from
        root (str): Directory to place downloaded file in
        filename (str, optional): Name to save the file under. If None, use the basename of the URL
    """
    result = 'OK'

    if not filename:
        filename = os.path.basename(url)
    fpath = os.path.join(root, filename)

    try:
        print('Downloading ' + url + ' to ' + fpath)
        request.urlretrieve(url, fpath)

    except (error.HTTPError, error.URLError, IOError) as e:

# Monitoring error 404 (file not found on rcsb website)
        if e.code == 404:
            print('=================> Erreur HTTP' +
                  ' while downloading ' + url + '  :  ' + f"{e.code}" + ' ' + e.reason)
            return ("ERR_404")

# Monitoring error 400 (file not found with https, trying with http)
        if e.code == 400:
            if url[:5] == 'https':
                url = url.replace('https:', 'http:')
                print('Failed download. Trying https -> http instead.'
                      ' Downloading ' + url + ' to ' + fpath)
                try:
                    request.urlretrieve(url, fpath)
                except (error.HTTPError, error.URLError, IOError) as e:
                    if e.code == 404:
                        print('=================> Erreur HTTP' +
                              ' while downloading ' + url + '  :  ' + f"{e.code}" + ' ' + e.reason)
                finally:
                    print("passage au fichier suivant")
    finally:
        print("passage au fichier suivant")

# Defining the function to read a text file
def read_text_file(file_path):
    fastaMatrix.clear()
    with open(file_path, 'r') as fn:
        i = 0
        for line in fn:

            if file_path.endswith("fasta"):
                if line.startswith('>'):
                    x = line.find('auth') + 5
                    if x == 4:
                        x = line.find('Chains') + 7
                    if x == 6:
                        x = line.find('Chain') + 6
# Initialization of the matrix containing in first column the chain nams and in second column the chains
                    fastaMatrix.append([line[line.index('|') + 1: line.index('|', line.index('|') + 1)].replace(
                        "Chains", "").replace("Chain", "").strip()])
                    header.write(f'{fasta_file} : {fastaMatrix[i][0]} \n')

                else:
                    fastaMatrix[i].append(line[0:].strip())
                    i += 1

# Defining the generator of matrix of common chains beetween the Fasta matrix and pdb matrix
def chercheChaineCommune(C1, C2):
    d1 = 0
    f1 = 0
    d2 = 0
    f2 = 0
    l = 2  # Common string = minimum 2 characters
    last_pos = 0

    while ((f1 < len(C1)) & (f2 < len(C2))):
        while (C1[f1:f1 + l] != C2[f2:f2 + l]) & (f1 < len(C1)) & (f2 < len(C2)):
            f1 += 1
        d1 = f1
        d2 = f2

        while ((C1[d1:d1 + l] == C2[d2:d2 + l]) & (d1 + l <= len(C1)) & (d2 + l <= len(C2))):
            l += 1
        if (C1[d1:d1 + l - 1] != ''):
            yield [d1, d1 + l - 2, C1[d1:d1 + l - 1]]
            last_pos = d1 + l - 1

        if ((f1 >= len(C1) - 1) & (f2 < len(C2) - 1)):
            f1 = last_pos
            f2 += 1
        else:
            f1 = d1 + l - 1
            f2 = d2 + l - 1
        l = 2


# ===========================================================================

# Main

# Variable declaration

Ligne = ""
result = ""
commonSubstrings = []

# Empty directory of downloaded FASTA files
path = "/home/erwan/Documents/Projet_Long/Fasta"
os.chdir(path)
for file in os.listdir():
    os.remove(file)

# Empty the directory of PDB files translated into Fasta
path = "/home/erwan/Documents/Projet_Long/Pdb2Fasta"
os.chdir(path)
for file in os.listdir():
    os.remove(file)

# Replacing the csv file
file_path = ''
if os.path.exists("/home/erwan/Documents/Projet_Long/resultats.csv"):
    os.remove("/home/erwan/Documents/Projet_Long/resultats.csv")

csv = open("/home/erwan/Documents/Projet_Long/resultats.csv", "a")
csv.write(f"Protéine;Chaine;Trous\n")

header = open("/home/erwan/Documents/Projet_Long/controle.txt", "w")

# Scanning the PDB directory
path = "/home/erwan/Documents/Projet_Long/NR_All_Martin"
os.chdir(path)

# iterate through all file

for file in os.listdir():
    # Check whether file is in pdb format or not
    if file.endswith(".pdb"):
        print("==============================================================================================")
        print("PDBFile  " + file)
        print("==============================================================================================")

        # Downloading corresponding Fasta File
        file_path = f"{path}/{file}"

        url = "https://www.rcsb.org/fasta/entry/" + file[:4] + "/"
        target = "/home/erwan/Documents/Projet_Long/Fasta"  # A rendre paramétrable
        fasta_file = file[:4] + ".fasta"
        result = download_url(url, target, fasta_file)

        if result == "ERR_404":
            Ligne = file[:4] + '; File not found error;\n'
            csv.write(Ligne)
            continue
        pdbMatrix.clear()

        PDBFile = file_path
        i = 0

        # Converting pdbFile into Fasta Format and loading in a matrix containing in the first column the identifier
        # and in the second the sequences
        pdb2fasta = open("/home/erwan/Documents/Projet_Long/Pdb2Fasta/" + fasta_file, "a")

        with open(PDBFile, 'r') as pdb_file:
            for record in SeqIO.parse(pdb_file, 'pdb-atom'):
                pdbMatrix.append([record.id[5:6]])
                pdb2fasta.write(">" + str([record.id[5:6]])[2:3] + '\n')
                ligne = record.seq
                pdb2fasta.write(str(ligne) + '\n')
                pdbMatrix[i].append(ligne)
                i += 1

        pdbMatrix.sort()

        # Loading fastaFile in memory
        read_text_file(f"{target}/{fasta_file}")

        fastaMatrix.sort()

        i = 0
        j = 0
        while i < len(pdbMatrix) and j < len(fastaMatrix):

            Ligne = file[:4] + ";"

            if len(pdbMatrix) < len(fastaMatrix):
                f = j
                p = -1
                for x in pdbMatrix:
                    p += 1
                    if x[0] in fastaMatrix[j][0]:
                        break

                if pdbMatrix[p][0] in fastaMatrix[j][0]:
                    i = p
                else:
                    Ligne += fastaMatrix[j][0] + ";Corresponding chain not found"
                    csv.write(f"{Ligne}\n")
                    j = j + 1
                    continue

     #       print('==========')
     #       print(f'fastaMatrix({j})')
     #       print(f'> {fastaMatrix[j][0]}')
     #       print(fastaMatrix[j][1])
     #       print(f'pdbMatrix({i})')
     #       print(f'> {pdbMatrix[i][0]}')
     #       print(pdbMatrix[i][1])
     #       print('============')

            if result == "ERR_404":
                Ligne = Ligne + "ERR_404;"
            else:
                pdb = pdbMatrix[i][1]
                fasta = fastaMatrix[j][1]

                Ligne = Ligne + pdbMatrix[i][0] + ";"

                commonSubstrings = []
                f = 0
                p = 0
                g = 0
                q = 0

                # Search for common subchains
                for substring in chercheChaineCommune(fasta, pdb):
                    commonSubstrings.append(substring)
                # print ("================== Sous-chaines Communes")
                # print(commonSubstrings)

                if len(commonSubstrings) == 0:
                    Ligne = file[:4] + '; No common chain or strings found;\n'
                    csv.write(f"{Ligne}\n")
                    i += 1
                    j += 1
                    continue

                c = 0
                if commonSubstrings[c][0] > 0:
                    Ligne += '0-' + str(commonSubstrings[c][0] - 1) + ';'

                c += 1

                while c < len(commonSubstrings):
                    Ligne += str(commonSubstrings[c - 1][1] + 1) + '-' + str(commonSubstrings[c][0] - 1) + ';'
                    c += 1

                if commonSubstrings[c - 1][1] < len(fasta) - 1:
                    Ligne += str(commonSubstrings[c - 1][1] + 1) + '-' + str(len(fasta) - 1)

                for x in commonSubstrings:
                    print(x)

                print(Ligne)

                csv.write(f"{Ligne}\n")
                i += 1
                j += 1
csv.close()
# RNK.read_data_set("/home/erwan/Documents/Projet_Long/resultats.csv")