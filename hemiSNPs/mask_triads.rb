#!/usr/bin/env ruby
require 'optparse'

require 'csv'
require 'fileutils'
require 'tmpdir'
require 'bio-samtools'
require 'bio'

class Array
  def sum
    inject(0.0) { |result, el| result + el }
  end

  def mean
    sum / size
  end
end

module Bio::PolyploidTools::Masks

  def self.find_end(seqs)
    size = seqs.values[0].size
    names = seqs.keys
    i = size - 1
    gap_count = 3
    while i > 0 and gap_count > 0
      gap_count = names.map { |chr| seqs[chr][i] == "-" ? 1:0  }.inject(0, :+)
      i -= 1
    end
    i + 1
  end

  def self.find_start(seqs)
    size = seqs.values[0].size
    names = seqs.keys
    i = 0
    gap_count = 3
    while i < size  and gap_count > 0
      gap_count = names.map { |chr| seqs[chr][i] == "-" ? 1 : 0  } .inject(0, :+)

      i += 1
    end
    i - 1
  end

  def self.get(seqs, target: nil, seq_start: 0, seq_end: 0)
    names = seqs.keys
    target = names[0] if target.nil?
    masked_snps = seqs[target].downcase
    i = 0
    while i < masked_snps.size
      different = 0
      cov = 0
      gap = false
      names.each do | chr |
        if seqs[chr][i]  != "-" and seqs[chr][i]  != "n" and seqs[chr][i]  != "N"
          cov += 1
        end
        if chr != target
         different += 1  if masked_snps[i].upcase != seqs[chr][i].upcase
        end
        if seqs[chr][i]  == "-" and chr == target
            gap = true
        end
      end
      masked_snps[i] = "." if different == 0
      masked_snps[i] = "." if cov == 1
      masked_snps[i] = "*" if cov == 0
      expected_snps  = names.size - 1
      masked_snps[i] = masked_snps[i].upcase if different == expected_snps
      if gap
        masked_snps[i] = different == expected_snps ? "-" : "_"
      end
      masked_snps[i] = "|" if i < seq_start or i > seq_end
      i += 1
    end
    masked_snps
  end

  def self.stats(mask, triad, gene, genome, reference)
    specific = []
    semispecific = []
    sp_i = 0
    semi = 0
    i = 0
    mask.to_s.each_char do |e|
      case e
      when "n","N"
        i += 1
      when /[[:lower:]]/ then
        semispecific << semi
        semi = 0
        i += 1
      when /[[:upper:]]/ then
        specific     << sp_i
        semispecific << semi
        sp_i = 0
        semi = 0
        i += 1
      when "." then
        semi += 1
        sp_i += 1
        i += 1
      end
    end
    {
      reference: reference,
      triad: triad,
      genome: genome,
      gene: gene,
      semispecific_mean: semispecific.mean,
      semispecific_bases: semispecific.size,
      semispecific_identity: (1 - (semispecific.size.to_f / i)) * 100 ,
      specific_mean: specific.mean,
      specific_bases: specific.size,
      specific_identity: (1 - (specific.size.to_f / i )) * 100,
      aligned_length: i,
      specific: specific,
      semispecific: semispecific
    }
  end
end

opts = {}
opts[:identity] = 50
opts[:min_bases] = 200
opts[:split_token] = "."
opts[:tmp_folder]  = Dir.mktmpdir
opts[:random_sample] = 0
opts[:output_folder] = "."

OptionParser.new do |o|

  o.banner = "Usage: mask_triads.rb [options]"

  o.on("-t", "--triads FILE", "CSV file with the gene triad names in the named columns 'A','B' and 'D' ") do |o|
    opts[:triads] = o
  end

  o.on("-f", "--fasta FILE" , "FASTA file containing all the possible sequences. ") do |o|
    opts[:fasta] = o
  end

  o.on("-s", "--split_token CHAR", "Character used to split the sequence name. The name will be evarything before this token on the name of the sequences") do |o|
    opts[:split_token] = o
  end

  o.on("-o", "--output_folder DIR", "Location to save the alignment masks. If the alignment exists, it is recycled to avoid calling MAFFT again") do |o|
    opts[:output_folder] = o
  end
end.parse!


split_token = opts[:split_token]
reference_name = File.basename opts[:fasta]
output_folder = opts[:output_folder]
@fasta_reference_db = Bio::DB::Fasta::FastaFile.new(fasta: opts[:fasta])
@fasta_reference_db.load_fai_entries
#puts @fasta_reference_db.index.entries
@cannonical = Hash.new
@fasta_reference_db.index.entries.each do |e|
  gene = e.id.split(split_token)[0]
  @cannonical[gene] = e unless @cannonical[gene]
  @cannonical[gene]  = e if   e.length > @cannonical[gene].length
end

$stderr.puts "#Loaded #{@cannonical.length} canonical sequences from #{@fasta_reference_db.index.size} in reference"

$stderr.puts "TMP dir: #{opts[:tmp_folder]}"

def write_fasta_from_hash(sequences, filename)
  out = File.new(filename, "w")
  sequences.each_pair do | chromosome, exon_seq |
    out.puts ">#{chromosome}\n#{exon_seq}\n"
  end
  out.close
end

def mafft_align(a, b, d)
  to_align = Bio::Alignment::SequenceHash.new
  seq_a = @fasta_reference_db.fetch_sequence(@cannonical[a].get_full_region)
  seq_b = @fasta_reference_db.fetch_sequence(@cannonical[b].get_full_region)
  seq_d = @fasta_reference_db.fetch_sequence(@cannonical[d].get_full_region)
  to_align[a] = seq_a
  to_align[b] = seq_b
  to_align[d] = seq_d
  report = mafft.query_alignment(to_align)
  aln = report.alignment
  aln
end

def read_alignment(path)
  aln = Bio::Alignment::SequenceHash.new
  i = 0
  Bio::FlatFile.open(Bio::FastaFormat, path) do |fasta_file|
    fasta_file.each do |entry|
      aln[entry.entry_id] = entry.seq if i < 3
      i += 1
    end
  end
  aln
end


mafft_opts = ['--maxiterate', '1000', '--localpair', '--quiet']
mafft = Bio::MAFFT.new( "mafft" , mafft_opts)
header_printed = false
stats     = File.open("#{output_folder}/#{reference_name}.identity_stats.csv", "w")
distances = File.open("#{output_folder}/#{reference_name}.distance_between_snps.csv.gz", "w")
gz = Zlib::GzipWriter.new(distances)
gz.write "triad,gene,genome,reference,type,distance\n"
#gz.close

def write_distances(distances, triad, gene, genome, reference, type, out)
  distances.each { |e| out.write "#{triad},#{gene},#{genome},#{reference},#{type},#{e}\n" }
end

i = 0
CSV.foreach(opts[:triads], headers:true ) do |row|
  next unless row["cardinality_abs"] == "1:1:1" and row["HC.LC"] == "HC-only"
   a = row['A']
   b = row['B']
   d = row['D']
   triad = row['group_id']
   cent_triad = triad.to_i / 100
   folder = "#{output_folder}/alignments/#{reference_name}/#{cent_triad}/"
   save_cds = "#{folder}/#{triad}.fa"
   aligned = File.file?(save_cds)
   aln = aligned ? read_alignment(save_cds)  : mafft_align(a,b,d)
   folder = "#{output_folder}/alignments_new/#{reference_name}/#{cent_triad}/" if aligned
   FileUtils.mkdir_p folder
   save_cds = "#{folder}/#{triad}.fa"

   aln2 = Bio::Alignment.new aln
   seq_start = Bio::PolyploidTools::Mask.find_start(aln)
   seq_end   = Bio::PolyploidTools::Mask.find_end(aln)
   #puts "#{triad}: #{seq_start}-#{seq_end}"


   aln2.add_seq(Bio::PolyploidTools::Mask.get(aln,seq_start: seq_start, seq_end: seq_end, target: a), "A")
   aln2.add_seq(Bio::PolyploidTools::Mask.get(aln,seq_start: seq_start, seq_end: seq_end, target: b), "B")
   aln2.add_seq(Bio::PolyploidTools::Mask.get(aln,seq_start: seq_start, seq_end: seq_end, target: d), "D")

   a_stats =  Bio::PolyploidTools::Mask.stats(aln2["A"], triad, a, "A", reference_name)
   b_stats =  Bio::PolyploidTools::Mask.stats(aln2["B"], triad, b, "B", reference_name)
   d_stats =  Bio::PolyploidTools::Mask.stats(aln2["D"], triad, d, "D", reference_name)
  
   write_distances(a_stats[:specific], triad, a, "A", reference_name, "specific", gz)
   write_distances(b_stats[:specific], triad, b, "B", reference_name, "specific", gz)
   write_distances(d_stats[:specific], triad, d, "D", reference_name, "specific", gz)

   write_distances(a_stats[:semispecific], triad, a, "A", reference_name, "semispecific", gz)
   write_distances(b_stats[:semispecific], triad, b, "B", reference_name, "semispecific", gz)
   write_distances(d_stats[:semispecific], triad, d, "D", reference_name, "semispecific", gz)
   
   a_stats.delete(:semispecific)
   b_stats.delete(:semispecific)
   d_stats.delete(:semispecific)
   
   a_stats.delete(:specific)
   b_stats.delete(:specific)
   d_stats.delete(:specific)

   a_stats[:length] = @cannonical[a].length
   b_stats[:length] = @cannonical[b].length
   d_stats[:length] = @cannonical[d].length

   stats.puts a_stats.keys.join(",") unless header_printed
   stats.puts a_stats.values.join(",")
   stats.puts b_stats.values.join(",")
   stats.puts d_stats.values.join(",")
   header_printed = true

   write_fasta_from_hash(aln2, save_cds)
   i += 1
end
gz.close
distances.close
stats.close
