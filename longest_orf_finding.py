import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
def analyze_sequence(record, file):
    seq= record.seq.upper()
    content_A=seq.count('A')
    content_G=seq.count('G')
    content_C=seq.count('C')
    content_T=seq.count('T')
    g_c_content= (content_G + content_C) / len(seq) * 100
    return{'file_stem': file.stem,
           'id': record.id,
          'gc_content': g_c_content,
          'length': len(seq),
          'description':record.description,
          'A':content_A,
          'G':content_G,
          'C':content_C,
          'T':content_T}
def find_orf_nucleotides(sequence):
    orf_list=[]
    stop_codons = ['TAA', 'TGA', 'TAG']
    seq_str = str(sequence)
    start_pos=0
    while True:
        start_pos = seq_str.find('ATG', start_pos)
        if start_pos == -1:
            break
        for stop_pos in range(start_pos + 3, len(seq_str) - 2, 3):
            codon = seq_str[stop_pos:stop_pos + 3]
            if codon in stop_codons:
                orf_list.append((start_pos, stop_pos + 3))
                break
        start_pos += 1
    return (orf_list)
def find_longest_orf(sequence):
    orf_list=find_orf_nucleotides(sequence)
    if not orf_list:
        return None,'', 0, '', 0, []
    longest_coordinates = None
    longest_dna = ''
    longest_dna_length = 0
    longest_protein = ''
    longest_protein_length = 0
    all_orfs=[]
    for start, end in orf_list:
        orf_dna = sequence[start:end]
        current_length = len(orf_dna)
        protein = str(Seq(orf_dna).transcribe().translate(to_stop=True))
        protein_len = len(protein)
        all_orfs.append({
            'start': start,
            'end': end,
            'dna_length': current_length,
            'protein_length': protein_len,
            'protein': protein + '...' if len(protein) > 50 else protein
        })
        if current_length > longest_dna_length:
                longest_dna_length= current_length
                longest_coordinates= (start, end)
                longest_dna = orf_dna
                longest_protein= protein
                longest_protein_length = protein_len
    return(longest_coordinates,longest_dna, longest_dna_length, longest_protein, longest_protein_length, all_orfs)
def find_protein_molecular_weight(protein_sequence):
        return(molecular_weight(protein_sequence, seq_type='protein'))
folder= Path('data')
#folder.glob('*.fasta')
results=[]
for file in folder.glob('*.fasta'):
    print(file.name)
    record= SeqIO.read(file, 'fasta')
    seq = record.seq.upper()
    coordinates, dna, dna_length, protein, protein_length,all_orfs = find_longest_orf(seq)
    
    if coordinates is not None:
        mw=find_protein_molecular_weight(protein)
        start, stop = coordinates
        print(f"  Coordinates: {start} to {stop} (nucleotide positions)")
        print(f"  DNA length: {dna_length} bp")
        print(f"  Protein length: {protein_length} amino acids")
        print(f"  Protein: {protein}" if len(protein) > 50 else f"  Protein: {protein}")

    else:
        mw= None
        print("No ORFs found!")
    result= analyze_sequence(record, file)
    result["longest_coordinates"] = coordinates
    result["longest_dna"] = str(dna)[:100] + "..." if dna and len(dna) > 100 else str(dna)
    result["longest_dna_length"] = dna_length
    result["longest_protein"] = protein + "..." if protein and len(protein) > 100 else protein
    result["longest_protein_length"] = protein_length
    result["total_orfs"] = len(all_orfs)
    result["all_orfs"] = all_orfs[:50]
    result['protein_MW']=mw
    results.append(result)
for result in results:
    print(f"file.stem: {result['file_stem']}")
    print("-" * 60)
    print(f"id: {result['id']}")
    print("-" * 60)
    print(f"description: {result['description']}")
    print("-" * 60)
    print(f"length: {result['length']}")
    print("-" * 60)
    print(f"A: {result['A']}")
    print("-" * 60)
    print(f"T: {result['T']}")
    print("-" * 60)
    print(f"G: {result['G']}")
    print("-" * 60)
    print(f"C: {result['C']}")
    print("-" * 60)
    print(f"gc: {result['gc_content']:.2f}%")
    print("-" * 60)
    print(f"Longest coordinates: {result['longest_coordinates']}")
    print("-" * 60)
    print(f"Longest DNA length: {result['longest_dna_length']} bp")
    print("-" * 60)
    print(f"Longest Protein length: {result['longest_protein_length']} amino acids")
    print("-" * 60)
    print(f"Longest Protein: {result['longest_protein']}")
    print("-" * 60)
    print(f"protein molecular weight: {result['protein_MW']}")
if results:
    with open('project_sequence_analysis.csv','w', newline='') as csv_file:
        writer= csv.DictWriter(csv_file, fieldnames=results[0].keys())
        writer.writeheader()
        for result in results:
#--------------------convert lists n tuples to strings for CSV file-----------
            for key, value in result.items():
                if isinstance(value, (list, tuple)):
                    result[key] = str(value)
            writer.writerow(result)
print('csv file has created successfully')
