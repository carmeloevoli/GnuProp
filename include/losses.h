#ifndef BENIAMINO_LOSSES_H
#define BENIAMINO_LOSSES_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

template <typename T> class TableLoader {
public:
  TableLoader(const std::string &filename) : filename_(filename) {}

  bool loadTable(std::vector<std::vector<T>> &table) {
    std::ifstream file(filename_);
    if (!file.is_open()) {
      std::cerr << "Error: Failed to open file " << filename_ << std::endl;
      return false;
    }

    std::string line;
    while (std::getline(file, line)) {
      std::vector<T> row;
      std::istringstream iss(line);
      std::string cell;

      while (std::getline(iss, cell, '\t')) {
        try {
          T value = convertFromString<T>(cell);
          row.push_back(value);
        } catch (const std::exception &e) {
          std::cerr << "Error: Invalid value in file " << filename_
                    << std::endl;
          return false;
        }
      }

      if (row.size() != 4) {
        std::cerr << "Error: Each row in the table must have 4 columns in file "
                  << filename_ << std::endl;
        return false;
      }

      table.push_back(row);
    }

    file.close();
    return true;
  }

private:
  std::string filename_;

  template <typename U> U convertFromString(const std::string &str) {
    std::istringstream iss(str);
    U value;
    iss >> value;
    if (!iss.eof() || iss.fail()) {
      throw std::invalid_argument(
          "Failed to convert string to the specified type");
    }
    return value;
  }
};

class Losses() {
public:
  Losses() {}
  virtual ~Losses() = default;

private:
  std::string m_filename = "data/SimProp_proton_losses.txt";
}

#endif