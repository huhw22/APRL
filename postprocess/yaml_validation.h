#ifndef FEL_POSTPROCESS_YAML_VALIDATION_H
#define FEL_POSTPROCESS_YAML_VALIDATION_H

#include <initializer_list>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

#include "yaml-cpp/yaml.h"

namespace postprocess_common
{
  inline std::runtime_error cardError(const YAML::Node& node,
                                      const std::string& message)
  {
    std::ostringstream result;
    if (node.Mark().line >= 0)
      result << "YAML line " << (node.Mark().line + 1) << ": ";
    result << message;
    return std::runtime_error(result.str());
  }

  inline void validateMapKeys(const YAML::Node& node,
                              const std::string& context,
                              std::initializer_list<const char*> keys)
  {
    if (!node || !node.IsMap())
      throw cardError(node, context + " must be a map");
    std::set<std::string> allowed;
    std::ostringstream choices;
    bool first = true;
    for (std::initializer_list<const char*>::const_iterator key = keys.begin();
         key != keys.end(); ++key)
      {
        allowed.insert(*key);
        if (!first) choices << ", ";
        choices << *key;
        first = false;
      }
    for (YAML::const_iterator entry = node.begin(); entry != node.end();
         ++entry)
      {
        std::string key;
        try { key = entry->first.as<std::string>(); }
        catch (const YAML::Exception&) {
          throw cardError(entry->first,
            context + " keys must be scalar strings");
        }
        if (allowed.find(key) == allowed.end())
          throw cardError(entry->first, context + " contains unknown key '" +
            key + "'; recognized keys are: " + choices.str());
      }
  }

  inline void validateSequenceMaps(const YAML::Node& sequence,
                                   const std::string& context,
                                   std::initializer_list<const char*> keys)
  {
    if (!sequence || !sequence.IsSequence())
      throw cardError(sequence, context + " must be a sequence");
    for (std::size_t i = 0; i < sequence.size(); ++i)
      {
        std::ostringstream item;
        item << context << '[' << i << ']';
        validateMapKeys(sequence[i], item.str(), keys);
      }
  }

  inline void requireOutputAvailable(const std::string& path, bool overwrite)
  {
    struct stat status;
    if (!overwrite && ::stat(path.c_str(), &status) == 0)
      throw std::runtime_error("Refusing to overwrite existing output: " +
        path + "; set output.overwrite: true only when replacement is intentional");
  }
}

#endif
