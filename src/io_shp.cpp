// [[Rcpp::depends(RcppParallel)]]

#include "io_shp.h"
#include <RcppParallel.h>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstdint>
using namespace RcppParallel;

namespace {

std::vector<unsigned char> read_whole_file(const std::string &path){
  std::ifstream f(path.c_str(), std::ios::binary | std::ios::ate);
  if (!f) Rcpp::stop("Cannot open file: " + path);
  std::streamsize size = f.tellg();
  f.seekg(0, std::ios::beg);
  std::vector<unsigned char> buf(static_cast<size_t>(size));
  if (size > 0 && !f.read(reinterpret_cast<char*>(buf.data()), size))
    Rcpp::stop("Failed reading file: " + path);
  return buf;
}

bool host_is_little_endian(){
  uint16_t x = 1;
  return *reinterpret_cast<unsigned char*>(&x) == 1;
}

int32_t read_be_i32(const unsigned char *p){
  return (static_cast<int32_t>(p[0]) << 24) | (static_cast<int32_t>(p[1]) << 16) |
         (static_cast<int32_t>(p[2]) << 8)  |  static_cast<int32_t>(p[3]);
}

int32_t read_le_i32(const unsigned char *p){
  static const bool le = host_is_little_endian();
  int32_t v;
  std::memcpy(&v, p, 4);
  if (!le) v = (static_cast<int32_t>(p[3]) << 24) | (static_cast<int32_t>(p[2]) << 16) |
              (static_cast<int32_t>(p[1]) << 8)  |  static_cast<int32_t>(p[0]);
  return v;
}

double read_le_f64(const unsigned char *p){
  static const bool le = host_is_little_endian();
  double v;
  if (le){
    std::memcpy(&v, p, 8);
  } else {
    unsigned char tmp[8];
    for (int i = 0; i < 8; i++) tmp[i] = p[7 - i];
    std::memcpy(&v, tmp, 8);
  }
  return v;
}

// Fixed layout offsets, relative to the start of a record's shape content
// (i.e. right after the 8-byte record header), for PolyLine/PolyLineZ/PolyLineM:
// +0 shape type (int32), +4 box (4 doubles = 32 bytes), +36 NumParts (int32),
// +40 NumPoints (int32), +44 Parts array (NumParts x int32), then Points
// array (NumPoints x (X,Y) doubles). Z/M arrays (if any) follow the points
// and are skipped entirely -- graph topology only needs X/Y.
const size_t PARTS_OFFSET = 44;

struct ShpGeomWorker : public Worker {
  const std::vector<unsigned char> &buf;
  const std::vector<size_t> &rec_content_start;
  const std::vector<int> &rec_num_parts;
  const std::vector<size_t> &part_row_offset;

  std::vector<double> &out_x1;
  std::vector<double> &out_y1;
  std::vector<double> &out_x2;
  std::vector<double> &out_y2;
  std::vector<double> &out_length;
  std::vector<int> &out_record;

  ShpGeomWorker(const std::vector<unsigned char> &buf_,
                const std::vector<size_t> &rec_content_start_,
                const std::vector<int> &rec_num_parts_,
                const std::vector<size_t> &part_row_offset_,
                std::vector<double> &out_x1_, std::vector<double> &out_y1_,
                std::vector<double> &out_x2_, std::vector<double> &out_y2_,
                std::vector<double> &out_length_, std::vector<int> &out_record_) :
    buf(buf_), rec_content_start(rec_content_start_), rec_num_parts(rec_num_parts_),
    part_row_offset(part_row_offset_), out_x1(out_x1_), out_y1(out_y1_),
    out_x2(out_x2_), out_y2(out_y2_), out_length(out_length_), out_record(out_record_) {}

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t r = begin; r < end; r++){
      int num_parts = rec_num_parts[r];
      if (num_parts == 0) continue;

      const unsigned char *content = buf.data() + rec_content_start[r];
      int num_points = read_le_i32(content + 40);
      const unsigned char *parts = content + PARTS_OFFSET;
      const unsigned char *points = parts + 4 * static_cast<size_t>(num_parts);

      for (int k = 0; k < num_parts; k++){
        int part_start = read_le_i32(parts + 4 * k);
        int part_end = (k + 1 < num_parts) ? read_le_i32(parts + 4 * (k + 1)) : num_points;
        if (part_end <= part_start) continue;

        size_t row = part_row_offset[r] + static_cast<size_t>(k);

        double x1 = read_le_f64(points + 16 * static_cast<size_t>(part_start));
        double y1 = read_le_f64(points + 16 * static_cast<size_t>(part_start) + 8);
        double x2 = read_le_f64(points + 16 * static_cast<size_t>(part_end - 1));
        double y2 = read_le_f64(points + 16 * static_cast<size_t>(part_end - 1) + 8);

        double length = 0.0;
        double px = x1, py = y1;
        for (int j = part_start + 1; j < part_end; j++){
          double qx = read_le_f64(points + 16 * static_cast<size_t>(j));
          double qy = read_le_f64(points + 16 * static_cast<size_t>(j) + 8);
          length += std::sqrt((qx - px) * (qx - px) + (qy - py) * (qy - py));
          px = qx; py = qy;
        }

        out_x1[row] = x1; out_y1[row] = y1;
        out_x2[row] = x2; out_y2[row] = y2;
        out_length[row] = length;
        out_record[row] = static_cast<int>(r);
      }
    }
  }
};

std::string trim(const char *p, size_t len){
  size_t a = 0, b = len;
  while (a < b && (p[a] == ' ' || p[a] == '\0')) a++;
  while (b > a && (p[b-1] == ' ' || p[b-1] == '\0')) b--;
  return std::string(p + a, b - a);
}

// .dbf files predate UTF-8 and are commonly written in a single-byte Windows
// codepage -- CP1252 by far the most common case for European (esp. French
// government/IGN) shapefiles. R sessions are UTF-8 by default nowadays, so a
// raw byte passthrough silently breaks equality comparisons against any
// accented literal (e.g. attrs$ACCES == "A péage" never matches). Detect the
// codepage via the sidecar .cpg file (or, failing that, the dbf header's
// legacy "language driver" byte) and convert to UTF-8 when it says CP1252.
bool codepage_is_cp1252(const std::string &dbf_path, const std::vector<unsigned char> &dbf_buf){
  std::string cpg_path = dbf_path.substr(0, dbf_path.find_last_of('.')) + ".cpg";
  std::ifstream cpg(cpg_path.c_str());
  if (cpg){
    std::string tag; std::getline(cpg, tag);
    return tag.find("1252") != std::string::npos;
  }
  // Fallback: dbf header language driver byte (offset 29). 0x03 is the
  // long-standing convention for "Windows ANSI codepage" (CP1252).
  return dbf_buf.size() > 29 && dbf_buf[29] == 0x03;
}

void utf8_append(unsigned int cp, std::string &out){
  if (cp < 0x80){ out += static_cast<char>(cp); }
  else if (cp < 0x800){
    out += static_cast<char>(0xC0 | (cp >> 6));
    out += static_cast<char>(0x80 | (cp & 0x3F));
  } else {
    out += static_cast<char>(0xE0 | (cp >> 12));
    out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
    out += static_cast<char>(0x80 | (cp & 0x3F));
  }
}

std::string cp1252_to_utf8(const std::string &in){
  // WHATWG-standard CP1252 code points for the 0x80-0x9F range (undefined
  // slots map to themselves, same as every major CP1252 implementation).
  static const unsigned int hi[32] = {
    0x20AC,0x0081,0x201A,0x0192,0x201E,0x2026,0x2020,0x2021,
    0x02C6,0x2030,0x0160,0x2039,0x0152,0x008D,0x017D,0x008F,
    0x0090,0x2018,0x2019,0x201C,0x201D,0x2022,0x2013,0x2014,
    0x02DC,0x2122,0x0161,0x203A,0x0153,0x009D,0x017E,0x0178
  };
  std::string out;
  out.reserve(in.size());
  for (unsigned char c : in){
    if (c < 0x80) out += static_cast<char>(c);
    else if (c < 0xA0) utf8_append(hi[c - 0x80], out);
    else utf8_append(c, out); // 0xA0-0xFF: CP1252 code point == byte value (Latin-1 supplement)
  }
  return out;
}

}

Rcpp::List read_shp_geometry(std::string path){
  std::vector<unsigned char> buf = read_whole_file(path);
  if (buf.size() < 100) Rcpp::stop("Not a valid .shp file (too small): " + path);

  int32_t file_code = read_be_i32(buf.data());
  if (file_code != 9994) Rcpp::stop("Not a valid .shp file (bad file code): " + path);

  int32_t shape_type = read_le_i32(buf.data() + 32);
  if (shape_type != 3 && shape_type != 13 && shape_type != 23)
    Rcpp::stop("Unsupported shape type (only PolyLine/PolyLineZ/PolyLineM are supported): " +
               std::to_string(shape_type));

  std::vector<size_t> rec_content_start;
  std::vector<int> rec_num_parts;
  std::vector<size_t> part_row_offset;
  size_t total_parts = 0;

  size_t offset = 100;
  while (offset + 8 <= buf.size()){
    int32_t content_length_words = read_be_i32(buf.data() + offset + 4);
    size_t content_length_bytes = static_cast<size_t>(content_length_words) * 2;
    size_t content_start = offset + 8;
    if (content_start + content_length_bytes > buf.size()) break;

    int32_t rec_shape_type = read_le_i32(buf.data() + content_start);
    int num_parts = 0;
    if (rec_shape_type == 3 || rec_shape_type == 13 || rec_shape_type == 23){
      num_parts = read_le_i32(buf.data() + content_start + 36);
    }

    rec_content_start.push_back(content_start);
    rec_num_parts.push_back(num_parts);
    part_row_offset.push_back(total_parts);
    total_parts += static_cast<size_t>(num_parts);

    offset = content_start + content_length_bytes;
  }

  std::vector<double> out_x1(total_parts), out_y1(total_parts), out_x2(total_parts),
                       out_y2(total_parts), out_length(total_parts);
  std::vector<int> out_record(total_parts);

  ShpGeomWorker worker(buf, rec_content_start, rec_num_parts, part_row_offset,
                        out_x1, out_y1, out_x2, out_y2, out_length, out_record);
  parallelFor(0, rec_content_start.size(), worker);

  return Rcpp::List::create(
    Rcpp::Named("x1") = Rcpp::wrap(out_x1),
    Rcpp::Named("y1") = Rcpp::wrap(out_y1),
    Rcpp::Named("x2") = Rcpp::wrap(out_x2),
    Rcpp::Named("y2") = Rcpp::wrap(out_y2),
    Rcpp::Named("length") = Rcpp::wrap(out_length),
    Rcpp::Named("record") = Rcpp::wrap(out_record),
    Rcpp::Named("n_records") = static_cast<int>(rec_content_start.size())
  );
}

Rcpp::DataFrame read_dbf(std::string path){
  std::vector<unsigned char> buf = read_whole_file(path);
  if (buf.size() < 32) Rcpp::stop("Not a valid .dbf file (too small): " + path);

  int32_t num_records = read_le_i32(buf.data() + 4);
  int16_t header_size; std::memcpy(&header_size, buf.data() + 8, 2);
  int16_t record_size; std::memcpy(&record_size, buf.data() + 10, 2);

  struct Field { std::string name; char type; int length; size_t offset; };
  std::vector<Field> fields;
  size_t desc_offset = 32;
  size_t field_offset = 1; // skip per-record deletion flag byte
  while (desc_offset + 32 <= static_cast<size_t>(header_size) && buf[desc_offset] != 0x0D){
    std::string name = trim(reinterpret_cast<const char*>(buf.data() + desc_offset), 11);
    char type = static_cast<char>(buf[desc_offset + 11]);
    int length = buf[desc_offset + 16];
    fields.push_back(Field{name, type, length, field_offset});
    field_offset += static_cast<size_t>(length);
    desc_offset += 32;
  }

  bool is_cp1252 = codepage_is_cp1252(path, buf);

  Rcpp::List cols(fields.size());
  Rcpp::CharacterVector names(fields.size());

  for (size_t fi = 0; fi < fields.size(); fi++){
    const Field &fld = fields[fi];
    names[fi] = fld.name;
    bool numeric = (fld.type == 'N' || fld.type == 'F');

    if (numeric){
      Rcpp::NumericVector col(num_records);
      for (int r = 0; r < num_records; r++){
        size_t row_start = static_cast<size_t>(header_size) + static_cast<size_t>(r) * record_size;
        std::string txt = trim(reinterpret_cast<const char*>(buf.data() + row_start + fld.offset), fld.length);
        col[r] = txt.empty() ? NA_REAL : std::atof(txt.c_str());
      }
      cols[fi] = col;
    } else {
      Rcpp::CharacterVector col(num_records);
      for (int r = 0; r < num_records; r++){
        size_t row_start = static_cast<size_t>(header_size) + static_cast<size_t>(r) * record_size;
        std::string txt = trim(reinterpret_cast<const char*>(buf.data() + row_start + fld.offset), fld.length);
        col[r] = is_cp1252 ? cp1252_to_utf8(txt) : txt;
      }
      cols[fi] = col;
    }
  }

  cols.attr("names") = names;
  cols.attr("class") = "data.frame";
  cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -num_records);
  return cols;
}

// [[Rcpp::export]]
Rcpp::List cpp_read_shp_geometry(std::string path){
  return read_shp_geometry(path);
}

// [[Rcpp::export]]
Rcpp::DataFrame cpp_read_dbf(std::string path){
  return read_dbf(path);
}
