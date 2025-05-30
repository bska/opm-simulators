/*
  Copyright 2021 Equinor.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_KEYWORDVALIDATION_HEADER_INCLUDED
#define OPM_KEYWORDVALIDATION_HEADER_INCLUDED

#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/common/OpmLog/KeywordLocation.hpp>
#include <opm/simulators/flow/ValidationFunctions.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm
{

class Deck;
class DeckKeyword;
class ErrorGuard;
class ParseContext;

namespace KeywordValidation
{
    // Describe an unsupported keyword:
    struct UnsupportedKeywordProperties {
        bool critical; // Set to true if presence of the keyword should be an error
        std::optional<std::string> message; // An optional message to show if the keyword is present
    };

    // Describe a partially or fully supported keyword item, by listing legal values:
    template <typename T>
    struct SupportedKeywordProperties {
        bool critical; // Set to true if an unsupported or invalid item value should be an error
        std::function<bool(T)> validator; // Predicate function to test values
        std::optional<std::string> message; // An optional message to show if an illegal item is encountered
    };

    // This is used to list unsupported kewyords.
    using UnsupportedKeywords = std::map<std::string, UnsupportedKeywordProperties>;

    // This is used to list the partially supported items of a keyword:
    template <typename T>
    using SupportedSingleKeywordItems = std::map<std::size_t, SupportedKeywordProperties<T>>;

    // This is used to list the keywords that have partially supported items or items that benefit from early validation:
    template <typename T>
    using SupportedKeywordItems = std::map<std::string, SupportedSingleKeywordItems<T>>;

    // This contains the information needed to report a single error occurence.
    // The validator will construct a vector of these, copying the relevant
    // information from the properties of the offending keywords and items.
    struct ValidationError {
        bool critical; // Determines if the encountered problem should be an error or a warning
        KeywordLocation location; // Location information (keyword name, file and line number)
        std::size_t record_number; // Number of the offending record
        std::optional<std::size_t> item_number; // Number of the offending item
        std::optional<std::string> item_value; // The offending value of a problematic item
        std::optional<std::string> user_message; // An optional message to show if a problem is encountered
    };

    // Get a formatted error report from a vector of validation errors. Set
    // include_noncritical to true if the report should include noncritical errors, and
    // include_critical to true if the report should include critical errors. These may
    // be set independently. If no errors are included the result will be an empty string.
    std::string get_error_report(const std::vector<ValidationError>& errors,
                                 const bool include_noncritical,
                                 const bool include_critical);

struct SupportedKeywords {
    const SupportedKeywordItems<std::string> string_items;
    const SupportedKeywordItems<int> int_items;
    const SupportedKeywordItems<double> double_items;
};

    class KeywordValidator
    {
    public:
        KeywordValidator(const UnsupportedKeywords& unsupported_keywords,
                         const SupportedKeywords& partially_supported_keywords,
                         const SupportedKeywords& fully_supported_keywords,
                         const std::unordered_map<std::string, ValidationFunction>& special_validation)
            : m_unsupported_keywords(unsupported_keywords)
            , m_partially_supported_keywords(partially_supported_keywords)
            , m_fully_supported_keywords(fully_supported_keywords)
            , m_special_validation(special_validation)
        {
        }

        // Validate a deck, reporting warnings and errors. If there are only
        // warnings, these will be reported. If there are errors, these are
        // reported, and execution of the program is halted, unless the argument
        // treat_critical_as_noncritical is true, then these also will only be
        // reported and not cause termination.
        void validateDeck(const Deck& deck,
                          const ParseContext& parse_context,
                          const bool treat_critical_as_noncritical,
                          ErrorGuard& error_guard) const;

        // Validate a single deck keyword. If a problem is encountered, add the
        // relevant information to the errors vector.
        void validateDeckKeyword(const DeckKeyword& keyword, std::vector<ValidationError>& errors) const;

    private:
        template <typename T>
        void validateKeywordItem(const DeckKeyword& keyword,
                                 const SupportedKeywordProperties<T>& properties,
                                 const bool multiple_records,
                                 const std::size_t record_number,
                                 const std::size_t item_number,
                                 const T& item_value,
                                 std::vector<ValidationError>& errors) const;

        void validateKeywordItems(const DeckKeyword& keyword,
                                  const SupportedKeywords& keyword_items,
                                  std::vector<ValidationError>& errors) const;

        template <typename T>
        void validateKeywordItems(const DeckKeyword& keyword,
                                  const SupportedKeywordItems<T>& supported_options,
                                  std::vector<ValidationError>& errors) const;

        const UnsupportedKeywords m_unsupported_keywords;
        const SupportedKeywords m_partially_supported_keywords;
        const SupportedKeywords m_fully_supported_keywords;
        const std::unordered_map<std::string, ValidationFunction> m_special_validation;
    };


    // Helper class to test if a given value is with a list of allowed values.
    template <typename T>
    class allow_values
    {
    public:
        allow_values(const std::initializer_list<T>& allowed_values)
        {
            std::copy(allowed_values.begin(),
                      allowed_values.end(),
                      std::back_inserter(m_allowed_values));
        }

        bool operator()(const T& value) const
        {
            return std::find(m_allowed_values.begin(), m_allowed_values.end(), value) != m_allowed_values.end();
        }

    private:
        std::vector<T> m_allowed_values;
    };

    // Helper to test if given string value is convertible to bool (see DeckItem::to_bool)
    struct is_bool_convertible {
        is_bool_convertible() {}
        bool operator()(const std::string& value) const {
            try {
                return DeckItem::to_bool(value) || true;
            } catch (const std::invalid_argument& e) {
                return false;
            }
        }
    };

} // namespace KeywordValidation

} // namespace Opm


#endif
